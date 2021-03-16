
![KVN](KvnFull.png)

# KVN - my pPXF helper 

pPXF Helper for use with BPASS - a test
This code is not open because pPXF is under copyright and also I can't be bothered to make it pretty

## Basic User Manual

### Making New Templates and saving them for posterity 
```
from KVN import KVN

kvn = KVN()

kvn.make_templates('../../BPASS_hoki_dev/bpass_v2.2.1_imf135_300/', fwhm_obs=FWHM_gal,  
                   wl_obs=wave, wl_range_obs=[WL[0],WL[-1]],
                   velscale_ratio=1, wl_range_padding=[-10,10],
                  ) 
                 
kvn.save('myKVN')           
```
| Parameter     | Description |
| ----------- | ----------- |
| `path_bpass_spectra` | Location of the folder where the bpass spectra are located|
| `wl_obs` | (Linear) Wavelength array of the observational data. If given, don't need log_wl_obs|
| `log_wl_obs` | (Natural Logarithm) Wavelength array of the observational data. If given, don't need wl_obs|
| `fwhm_obs` | FWHM of the observations. Applies if the dispersion is not wavelength dependent|
| `dispersion_obs` | Dispersion array. Same size ase wl_obs of log_wl_obs. Applies when the dispersion is wavelength dependent.|
| `wl_range_obs` | \[min, max\] Wavelength range of the observations|
| `velscale_ratio` |  Velocity scale ratio. Default is 1|
| `wl_range_padding` | \[-X,Y\] Padding given to the templates on creation. The template wl range will be \[wl_range_obs-X, wl_range_obs+Y\]. Default is \[-1,1\]|
| `binary` | Whether to include the binary model spectra. Default is True|
| `single` | Whether to include the single star model spectra. Default is False.|
| `z_list` | Which metallicities to consider, e.g. \['z020', 'z014'\]. Default is None. In which case all 13 BPASS metallicites are included.|
| `oversample` | The oversample keyword in ppxf |
| ` _max_age_index` | Maximum age index to include. Default is 42 which corresponds to a Hubble time.|

### pPXF solutions and KVN
```
pp = ppxf(kvn.templates, galaxy, noise, velscale, start,
          goodpixels=goodPixels, plot=True, moments=4,
          degree=4, vsyst=dv, velscale_ratio=kvn.velscale_ratio)
          
kvn.make_results(pp)
```

Then `kevin.results` will be a DataFramme with the ages metalicties and weights of the best fitting SEDs.

### Loading a pre-existing KVN object
```
from KVN import KVN

kvn = KVN()
kvn.load('myKVN')
```

## KVN Attributes 

### Defined by `make_templates()`

| Attributes | Description |
| ---------- | ------------|
| `.bpass_list_spectra` | List of BPASS spectral files |
| `.wl_obs` | (Linear) Wavelength array of the observational data. If given, don't need log_wl_obs|
| `.log_wl_obs` | (Natural Logarithm) Wavelength array of the observational data. If given, don't need wl_obs|
| `.wl_tem` | (Linear) Wavelength array of the BPASS templates. |
| `.log_wl_tem` | (Natural Logarithm) Wavelength array of the  BPASS templates. |
| `.dispersion_obs` | Dispersion array. Same size ase wl_obs of log_wl_obs. Applies when the dispersion is wavelength dependent.|
| `.velscale_obs` | Velocity scale of the observations (still a problem?) |
| `.velscale_tem` | Velocity scale of the templates. Same as `.velscale_obs` if `velscale_ratio` is set to 1 |
| `.wl_range_obs` | \[min, max\] Wavelength range of the observations |
| `.wl_range_tem` | \[min, max\] Wavelength range of the templates |
| `.velscale_ratio` |  Velocity scale ratio. Default is 1|
| `.fwhm_obs` | FWHM of the observations. Applies if the dispersion is not wavelength dependent|
| `.fwhm_tem` | FWHM of the templates |
| `.fwhm_dif`| Difference between obs and tem FWHM (in quadrature) |
| `.sigma` | `fwhm_dif/2.355/fwhm_tem` |
| `.z_list` | Which metallicities to consider, e.g. \['z020', 'z014'\]. Default is None. In which case all 13 BPASS metallicites are included.|
| `._max_age_index` | Maximum age index to include. Default is 42 which corresponds to a Hubble time. |
| `.templates` | **The templates** |
| `.template_preperties` | A matrix with dimensions (Nspectra, 2), where Nspectra = Nmetallicities * NAges. For each spectrum in the `templates` matrix, the metallcitiy and age is recorded in the `templates_properties` matrix |


### Defined by `make_results()`

| Attributes | Description |
| ---------- | ------------|
| `.ppxf` |  The ppxf object used to make results |
| `.matching_indices`| The indices of the matching templates and template parameters |
| `.matching_raw_spectra` | The raw BPASS spectra that match the observations | 
| `.matching_spectra` | The spectra combined to match observations convolved with the LOSVD from ppxf.matrix|
| `.matching_apolynomial` | The additive polynomial found by ppxf to help match observations. Can be None |
| `.matching_mpolynomial` | The multiplicative polynomial found by ppxf to help match observations. Can be None |
| `.results`| Data Frame containing the results: Age, Metallicity, Weights |

