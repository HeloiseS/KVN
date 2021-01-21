from hoki import load
from hoki.constants import BPASS_METALLICITIES, BPASS_NUM_METALLICITIES
import glob
from ppxf import ppxf_util
from os import path
import ppxf as ppxf_package
from astropy.io import fits
import numpy as np
from ppxf.ppxf import ppxf, losvd_rfft, rebin
from scipy import ndimage
import pandas as pd
import pickle
import matplotlib.pyplot as plt

from hoki.utils.hoki_object import HokiObject
from hoki.utils.hoki_dialogue import dialogue
from hoki.utils.exceptions import HokiFatalError
from hoki.utils.progressbar import print_progress_bar

plt.style.use('hfs')

####  Some useful definitions

met_to_num={'zem5':1e-5,  'zem4':1e-4, 'z001':1e-3, 'z002':2e-3,'z003':3e-3,'z004':4e-3,'z006':6e-3,'z008':8e-3,
            'z010':1e-2,'z014':1.4e-2,'z020':2e-2,'z030':3e-2, 'z040':4e-2}

c = 299792.458


class KVN(HokiObject):
    """
    Kevin - my pPXF helper
    We're not sure what he contributes but he's there for you
    """

    # TODO: Allow selection of ages??
    # isntanciated  (__ini__) by the hoki object

    def make_templates(self, path_bpass_spectra,
                       wl_obs=None, log_wl_obs=None,  fwhm_obs=None, dispersion_obs=None, wl_range_obs=None,
                       velscale_ratio=1, wl_range_padding=[-1,1],
                       binary=True, single = False, z_list=None, oversample=1, _max_age_index=42
                      ):

        print(f"{dialogue.info()} TemplateMaker Starting")
        print(f"{dialogue.running()} Initial Checks")

        self.model_path = path_bpass_spectra

        # Making list of relevant spectra files
        if binary and single:
            self.bpass_list_spectra = glob.glob(self.model_path + 'spectra*')
        elif binary and not single:
            self.bpass_list_spectra = glob.glob(self.model_path  + 'spectra*bin*')
        elif single and not binary:
            self.bpass_list_spectra = glob.glob(self.model_path  + 'spectra*sin*')
        else:
            raise HokiFatalError(f"binary and single set to False \n\n{dialogue.debugger()} "
                                 f"You must choose whether to include binary models, single models or both"
                                 f"\nAt least one of the following parameters must be True  when you instanciate:"
                                 f"\nbinary; single")

        self.bpass_list_spectra.sort()

        # Allow user to select a list of metallicities?
        if z_list is not None:
            self.z_list=z_list
            self.num_z_list=[met_to_num.get(key) for key in z_list]
        else:
            self.z_list=BPASS_METALLICITIES
            self.num_z_list=BPASS_NUM_METALLICITIES


        not_you=[]
        for filepath in self.bpass_list_spectra:
            if filepath[-8:-4] not in self.z_list:
                not_you.append(filepath)

        self.bpass_list_spectra = list(set(self.bpass_list_spectra)-set(not_you))


        #self.wl_range_tem = np.array((10**self.log_wl_obs)[[0,-1]].astype(int))
        if wl_obs is None and log_wl_obs is None:
            raise HokiFatalError(f"wavelength of observational data not provided\n\n{dialogue.debugger()}"
                                 f"At least one of the following parameters must be provided when you instanciate:"
                                 f"\nwl_obs (wavelength arr in linear space); "
                                 f"log_wl_obs (wavelength arr in log10 space)")
        elif wl_obs is not None and log_wl_obs is None:
            self.wl_obs = wl_obs
            self.log_wl_obs = np.log(self.wl_obs)
            self.wl_range_tem = np.array(self.wl_obs[[0,-1]].astype(int))

        elif log_wl_obs is not None and wl_obs is None:
            self.log_wl_obs = log_wl_obs
            self.wl_obs=np.e**self.log_wl_obs
            self.wl_range_tem = np.array(self.wl_obs[[0,-1]].astype(int))

        ### just to be sure we have enough wavelength coverage
        self.wl_range_tem+=wl_range_padding

        print(f"{dialogue.running()} Loading  model spectrum")
        # Load one spectrum
        _ssp=load.model_output(self.bpass_list_spectra[0])
        _ssp.index=_ssp.WL # turning  the WL column into an index so crop the dataframe easily with WL
        _ssp.drop('WL', axis=1, inplace=True)

        self.wl_tem = np.arange(self.wl_range_tem[0], self.wl_range_tem[-1]+1) # is this the same as ssp.WL?

        _ssp60=_ssp.loc[self.wl_range_tem[0]:self.wl_range_tem[-1]]['6.0'].values
        #one spectrum at log(age)=6.0: needed to get log_wl_template, velscale_template

        self.velscale_ratio=velscale_ratio

        self.fwhm_tem = 1 # known for BPASS, res is 1 Angstrom

        ## If disperion given
        if dispersion_obs is not None:
            print(f"{dialogue.running()} Calulating obs. velocity scale")
            # Calc. velocity scale
            _frac = self.wl_obs[1]/self.wl_obs[0]
            _dwl_obs = (_frac - 1)*self.wl_obs
            self.velscale_obs = np.log(_frac)*c #speed of light # Velocity scale in km/s per pixel

            print(f"{dialogue.running()} Calulating FWHM")
            # Calc. FWHM galaxy (res.)
            _fwhm_obs = 2.355*dispersion_obs*_dwl_obs
            self.fwhm_obs = np.interp(self.wl_tem, self.wl_obs, _fwhm_obs)
            self.fwhm_dif = np.sqrt((self.fwhm_obs**2 - self.fwhm_tem**2).clip(0))

        elif dispersion_obs is None and fwhm_obs is not None:
            self.wl_range_obs = wl_range_obs
            self.fwhm_obs = fwhm_obs
            print(f"{dialogue.running()} Calulating obs. velocity scale -- No dispersion")
            __, __, self.velscale_obs = ppxf_util.log_rebin(self.wl_range_obs, self.wl_obs)
            #self.velscale_obs=c*np.log(self.wl_obs[1]/self.wl_obs[0])
            self.fwhm_dif = np.sqrt((self.fwhm_obs**2 - self.fwhm_tem**2))

        elif dispersion_obs is None and fwhm_obs is None:
            raise HokiFatalError(f"dispersion_obs and fwhm_obs are both None \n\n{dialogue.debugger()}\
                    \nAt least one of the following parameters must be provided when you instanciate:\
                    \ndispersion_obs; fwhm_obs")

        print(f"{dialogue.running()} Calculating template wavelength (log rebin) and velocity scale")

        _ssp_temp, self.log_wl_tem, self.velscale_tem = ppxf_util.log_rebin(list(self.wl_range_tem),
                                                                            _ssp60, oversample=oversample,
                                                                            velscale=self.velscale_obs/self.velscale_ratio
                                                                            )

        self._max_age_index=_max_age_index
        print(f"{dialogue.running()} Calculating sigma")

        self.sigma = self.fwhm_dif/2.355/self.fwhm_tem

        print(f"{dialogue.running()} Instanciating container arrays")
        # instanciate empty arrays
        self.templates = np.empty((_ssp_temp.size, self._max_age_index, len(self.bpass_list_spectra)))
        # # WL bins, # Ages, # mets

        # TODO: find a better name
        self.template_properties =  np.zeros((len(self.bpass_list_spectra)*self._max_age_index, 2))

        print(f"{dialogue.running()} Compiling your templates")

        i = 0
        l = int(len(self.bpass_list_spectra))*self._max_age_index

        # index k (for the metallicity location in stars_tempaltes), and file path
        for k, path in enumerate(self.bpass_list_spectra):
            # string to identify the metallicity, like 'z040' above
            met = self.num_z_list[np.argwhere(np.array(self.z_list)==path[-8:-4])[0][0]]

            # load the file corresponding to that metallicity
            # getting the SSPs for metallciity at index k
            ssps_k=load.model_output(path)
            ssps_k.index=ssps_k.WL # make the wavelength the index, it makes things much easier
            ssps_k.drop('WL', axis=1, inplace=True) # drop the now unnecessary column
            cropped_ssps_k=ssps_k.loc[self.wl_range_tem[0]:self.wl_range_tem[-1]]

            # index j (age), and column '[age]' - only itterate till the current age of the universe
            for j, col in enumerate(cropped_ssps_k.columns[:self._max_age_index]):
                # get the simple star population for that age
                ssp_jk = cropped_ssps_k[col].values
                # gaussian filter to change the rsolution of the OG template (I think... check that)
                if isinstance(self.sigma, float):
                    ssp_jk = ndimage.gaussian_filter1d(ssp_jk, self.sigma)
                else:
                    ssp_jk = ppxf_util.gaussian_filter1d(ssp_jk, self.sigma) # not the scipy convolution because it can't do variable sigma

                # logarithm binning
                sspNew, __, __ = ppxf_util.log_rebin(self.wl_range_tem, ssp_jk,
                                                velscale=self.velscale_obs/self.velscale_ratio)
                # add to the template array
                self.templates[:, j, k] = sspNew/np.median(sspNew)
                self.template_properties[i, :] = [met, float(col)]
                i+=1

                print_progress_bar(i, l, prefix = 'Progress:', suffix = 'Complete', length = 50)

        print(f"{dialogue.complete()} Templates compiled successfully")

    def make_results(self, ppxf):
        self.ppxf=ppxf
        # indices of weights that are non-zero
        self.match_indices = np.argwhere(ppxf.weights != 0.0).T[0]

        # the non-zero weights
        weights = ppxf.weights[self.match_indices]
        self.matching_raw_spectra = []

        #iterating over every solution
        for i in range(self.template_properties[self.match_indices].shape[0]):
            # finding matching age
            i_age = np.argwhere(np.round(self.t,2)==self.template_properties[self.match_indices][i,1])[0]
            # finding matching metallicity
            i_met = np.argwhere(np.round(self.num_z_list,5)==self.template_properties[self.match_indices][i,0])[0]
            # compiling the spectra
            self.matching_raw_spectra.append(self.templates[:, i_age, i_met])

        self.matching_raw_spectra = np.array(self.matching_raw_spectra)
        self.matching_spectra = ppxf.matrix[:, ppxf.degree+1:].T[self.match_indices]
        self.matching_apolynomial = ppxf.apoly
        self.matching_mpolynomial = ppxf.mpoly

        # compiling summary of the results
        self.results = pd.DataFrame(np.vstack((self.template_properties[self.match_indices].T, weights.T)).T,
                                    columns=['met', 'age', 'weights'])

    def save(self, path):
        f = open(path, 'wb')
        pickle.dump(self.__dict__, f, 2)
        f.close()

    def load(self, path):
        f = open(path, 'rb')
        tmp_dict = pickle.load(f)
        f.close()

        self.__dict__.update(tmp_dict)