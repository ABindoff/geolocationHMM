# TwilightFree method of geolocation

A new method of geolocation from noisy light data using a hidden Markov model

For my Honours project I developed a new method of light-based geolocation that relies on the overall pattern of day and night rather than an explicit dependence on twilights as per the threshold, curve, or template-fit methods. For an introduction to geolocation
using recorded light, please see https://en.wikipedia.org/wiki/Light_level_geolocator<br><br>
This repository holds important data for my thesis and a paper which has been accepted detailing the 'twilight free' method.<br>
  
## Important Update: TwilightFree R package now available

This git is important as it provides a data archive including supplementary material for the published article. However, `TwilightFree` has been packaged for R and is available (with tutorials) from https://github.com/ABindoff/TwilightFree  

  
    
    

### Supplementary materials 

[twilight_free_algorith_example_ses.Rmd](https://github.com/ABindoff/geolocationHMM/blob/master/twilight%20free%20algorithm%20example%20ses.Rmd)  contains a script that replicates an example from the methods paper and calculates measures of accuracy and precision. 

[quick_start_guide.md](https://github.com/ABindoff/geolocationHMM/blob/master/quick_start_guide.md) provides a quick-start tutorial for users who want to use the (very small) set of functions to analyse their own light data from archival tags. *A tutorial will be added for users who wish to improve their location estimates by including SST data*

[parameter_sensitivity_new_data.md](https://github.com/ABindoff/geolocationHMM/blob/master/parameter_sensitivity_new_data.md) extensively tests the sensitivity of the shading/noise and movement parameters using data from a single animal, then uses the optimal parameters from that example on three more datasets, finally performing another PSA on an animal displaying markedly different behaviours to the first. *This analysis is useful because it is expected that users supply shading/noise and movement parameters from previous studies, preferably from double-tagging experiment. The analysis demonstrates that parameters determined from a double-tagging experiment can inform analyses of other tag data from the same species, but that optimal parameters are not necessary for useful results in biology and ecology.*


Known bug:

- where the `grid` (spatial map) crosses the equator it is possible for the MAP estimate to occur (somewhat obviously) on the wrong side of the equator on any given day (the algorithm picks a mathematically plausible but ecologically implausible solution). One quick fix is to break the data into parts and provide a different `grid` for each part. Alternatively, calling `essieRaster(fit_object)` and finding the max in the appropriate hemisphere also works. A solution is currently in development and testing.  

  
  [![DOI](https://zenodo.org/badge/61974427.svg)](https://zenodo.org/badge/latestdoi/61974427)


