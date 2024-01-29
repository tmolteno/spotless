# Spotless - Radio Interferometry Imaging Algorithm
    
This is a point-source deconvolution algorithm (part of the CLEAN family) that works without gridding. 
[http://www.iram.fr/IRAMFR/GILDAS/doc/html/map-html/node37.html]

This is essentially  a grid-free version of the Cotton-Schwab algorithm with a different convergence an optimization steps.
[Relaxing the isoplanatism assumption in self-calibration; applications to low-frequency radio interferometry]

It does not require W-projection and handles
non-coplanar antennas without difficulty. It also works on all-sky images just fine.

## How does Spotless work

The CLEAN algorithm is essentially deconvolution by repeated subtraction. I think this is silly, hence spotless. Spotless works by building up a model of the field of view in terms of point sources using model-fitting in visibility space. 

Spotless deconvolutes the measured visibilities $V(u,v,w)$ into a sum of $N$ point-source visibilities $V_P(\theta, \phi)$ where $\theta$ and $\phi$ are the co-ordinates of the point source: i.e.,

$$ V(u,v,w) = \sum_{i=0}^{N} A_i V_P(\theta_i, \phi_i) + V_r(u,v,w) $$

where the $A_i$ are the brightness of each point source, and $V_r$ are the residual visibilities.

Spotless has two algorithms for doing this. The first, like CLEAN, is sequential location of point sources:

```math
\begin{eqnarray*}
V(u,v,w) & = & A_0 V_P(\theta_0, \phi_0) + V_1(u,v,w) \\
V_1(u,v,w) & = & A_1 V_P(\theta_1, \phi_1) + V_2(u,v,w) \\
 ... \\
V_N(u,v,w) & = & A_N V_P(\theta_N, \phi_N) + V_{N+1}(u,v,w)
\end{eqnarray*}
```
At each step the new point source is located using a minimizer from the residuals at that step:
```math
\begin{eqnarray*}
  P_i & = & \min_{A, \theta, \phi} E(V_0) \\
      & = & \min_{A, \theta, \phi} E(V - A V_P(\theta, \phi))
\end{eqnarray*}
```
where $E(V)$ is the total energy in the visibilities. So we find the point source that minimizes the remaining energy.

The nice thing is that the energy can be calculated directly from the visibilities, and so no gridding is required at all, either in image space or in visibility space.

This is possible because the integral of the fourier transform (F.T.) of the visibilities can be calculated directly from the visibilities without the F.T. This is Parseval's Theorem, for the DFT it becomes
```math
\sum_{n=0}^{N-1} \left| x_n \right|^2 = \frac{1}{N} \sum_{k=0}^{N-1} \left| X_k \right|^2
```
Thus as the visibilities are the F.T of the sky brightness, the sum of the absolute value of the visibilities is proportional
to the energy in the 'image'
```math
E(V)  \propto  \sum_{k=0}^{N-1} v_i v_i^\star
```


### MultiSpotless 

Multispotless (--multimodel command line option) uses a better (but slower) sequential location. It builds up a multiple-point-source model as the algorithm progresses.


### Termination Criterion

Both spotless variants terminate when the power in the residual stops decreasing.

## Results

Dirty Image                |  Spotless Image
:-------------------------:|:-------------------------:
![](https://github.com/tmolteno/spotless/blob/main/img/gridless.png)  |  ![](https://github.com/tmolteno/spotless/blob/main/img/spotless.png)

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

## Changes

* 0.4.6 Clean up logging. Use --debug to switch it on. Remove deprecated calls to verbose in healpy
* 0.4.3 Use a common field of view parser with DiSkO
        Add --hdf <filename> option to save the output as an HDF
* 0.4.2 Read from measurement sets --ms
* 0.4.1 Use the disko sphere.
        Clean up unused code.
        Use harmonics from disko.
        Specify the --fov and --res as in disko
* 0.4.0 Move to github repository. Add to pypi. Use disko for utility functions.
        Add a --version CLI argument
* 0.3.0 Update to python3
* 0.3.3 Add a gridless binary to plot a GRIDLESS imaging, add a PDF export option
* 0.3.4 Fix bitrot in multispotless.
