# Diploma Project

```
Usage: VirtualRayLights.exe [ -s <scene_id> | -spp <samples_per_pixel> | -vrlps <vrl_per_sample> | -is <importance_sampling_strategy> | -b <bounces> | -sa <sigma_a> | -ss <sigma_s> | -g <henyey_greenstein> | -i <light_intensity> | -mr <medium_radiance> | -sr <surface_radiance> ]

        -s      selects a scene (default 0):
                        0       cornell box with directional lights
                        1       bidirectional VRLs
                        2       unidirectional VRLs
                        3       medium radiance test
                        4       surface radiance test
                        5       media cube
        -spp    samples per pixel
        -vrlps  VRLs per sample
        -is     selects the importance sampling techniques for the light / camera rays (default 4):
                        0       uniform / uniform
                        1       exponential / exponential
                        2       uniform / equi-angular
                        3       simple joint distribution
                        4       advanced joint distribution
                        5       product of Gaussian distributions
                        6       product of Cauchy distributions
                        7       Mixture Model with Cauchy distributions
                        8       Mixture Model with Gaussian distributions
                        9       Mixture Model with Hyperbolic Secant distributions
                        10      Mixture Model with Logistic distributions
        -b      number of allowed light bounces
        -sa     sigma_a
        -ss     sigma_s
        -g      henyey-greestein
        -i      the intensity of the light
        -mr     option for medium radiance (default true): true / false
        -sr     option for surface radiance (default true): true / false
```
