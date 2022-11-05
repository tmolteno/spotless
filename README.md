# Spotless - Radio Interferometry Imaging Algorithm
    
This is a point-source deconvolution algorithm (part of the CLEAN family) that works without gridding. 
[http://www.iram.fr/IRAMFR/GILDAS/doc/html/map-html/node37.html]

This is essentially  a grid-free version of the Cotton-Schwab algorithm with a different convergence an optimization steps.
[Relaxing the isoplanatism assumption in self-calibration; applications to low-frequency radio interferometry]

It does not require W-projection and handles
non-coplanar antennas without difficulty. It also works on all-sky images just fine.

Spotless works in a similar way to Hogbom's CLEAN algorithm, however 

 ![Dirty, Spotless and Residual][tart_image] 

For more information see the [TART Github repository](https://github.com/tmolteno/TART)

## Install Instructions

tart_tools is available from standard python package repositories. Try:

    pip install spotless
    
## Running it on live data

    spotless --api https://tart.elec.ac.nz/signal --display --show-sources
    gridless --api https://tart.elec.ac.nz/signal --display --show-sources

## TODO

* Add Gaussian Source Model
* Make explicit the antenna model (gain as a function of angular coordinates). We are assuming it is hemispherical here.
* Prove the relationship between power in the image and visibilty amplitudes. This might only work when the image 
  tends towards a random one. But this is OK since as we remove the sources the residual becomes more and more random.
* Run an MCMC on the multimodel option to estimate uncertainty in the model. Then use this uncertainty as a stopping criterion (when new model components no longer have certain amplitude or position)

## Author

* Tim Molteno (tim@elec.ac.nz)

## Development work
    
If you are developing this package, this should be installed using
```
	make develop
```
in which case changes to the source-code will be immediately available to projects using it.

    
[spotless_image]: https://github.com/tmolteno/TART/blob/master/doc/img/tart_image.jpg "TART All-Sky Image"


## Changes

* 0.4.0 Move to github repository. Add to pypi. Use disko for utility functions.
        Add a --version CLI argument
* 0.3.0 Update to python3
* 0.3.3 Add a gridless binary to plot a GRIDLESS imaging, add a PDF export option
* 0.3.4 Fix bitrot in multispotless.
