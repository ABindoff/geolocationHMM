# TwilightFree method of geolocation

A new method of geolocation from noisy light data using a hidden Markov model

As part of my Honours project I developed a new method of light-based geolocation that relies on the overall pattern of day and night,
rather than an explicit dependence on twilights as per the threshold, curve, or template-fit methods. For an introduction to geolocation
using recorded light, please see https://en.wikipedia.org/wiki/Light_level_geolocator<br><br>
This repository holds important data for my thesis and a paper which is to be submitted detailing the 'twilight free' method.<br>

Some important files, 

[twilight_free_algorith_example_ses.Rmd](https://github.com/ABindoff/geolocationHMM/blob/master/twilight%20free%20algorithm%20example%20ses.Rmd)  contains a script that replicates an example from the methods paper and calculates measures of accuracy and precision. 

[quick_start_guide.Rmd](https://github.com/ABindoff/geolocationHMM/blob/master/quick_start_guide.Rmd) provides a quick-start tutorial for users who want to use the (very small) set of functions to analyse their own light data from archival tags. *A tutorial will be added for users who wish to improve their location estimates by including SST data*

[parameter_sensitivity_new_data.Rmd](https://github.com/ABindoff/geolocationHMM/blob/master/parameter_sensitivity_new_data.Rmd) tests the sensitivity of the shading/noise and movement parameters using the example data from the methods paper, then uses the optimal parameters from that example on three more datasets, finally performing another PSA on an animal displaying markedly different behaviours to the first. *This analysis is useful because it is expected that users supply shading/noise and movement parameters from previous studies, preferably from double-tagging experiment. The analysis demonstrates that parameters determined from a double-tagging experiment can inform analyses of other tag data from the same species.*