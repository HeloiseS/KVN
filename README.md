
![KVN](KvnFull.png)

# KVN - my pPXF helper 

pPXF Helper for use with BPASS.
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
