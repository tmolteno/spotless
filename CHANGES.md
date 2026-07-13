# Changes

## 0.7.5

* Upgrade disko dependency to >=1.4.0 (many fixes since 1.1.1)
* Add --save-model-json FILE option to export point-source model
* Add LaTeX article, bibtex bibliography, and result images in doc/
* Fix MultiSpotless convergence: maxfev budget, retry on failure
* Document memory scaling: O(n_vis × n_pix) harmonic cache with table
* LaTeX documentation updated to match actual implementation

## 0.7.4

* Add --log FILE option to save deconvolution statistics
* Show deconvolution progress and final model summary
* Add optimizer stats (nfev, nit, converged) to step output

## 0.7.3

* Fix UnboundLocalError in CLI --version output
* Pre-warm disko harmonic cache (~8x speedup per step)
* Add performance tuning section to README
* Update copyright headers to 2022-2026
* Add comprehensive CLI usage docs to README

## 0.7.2

* Add LaTeX article in doc/ directory
* Switch inline bibliography to bibtex with spotless.bib

## 0.7.1

* Fix MultiSpotless optimizer (revert to Nelder-Mead adaptive)

## 0.7.0

* Migrate from Poetry to uv for package management and builds
* Fix test imports and deprecated API calls
* Fix add_source to use incremental subtraction
* Improve optimization: use peak amplitude as initial guess for optimizer
* Switch Spotless to L-BFGS-B (fast convergence on well-scaled 3-param problem)
* Keep MultiSpotless on Nelder-Mead with adaptive option (handles mixed-scale parameter space)

## 0.6.1

* Initial Poetry-based release

## 0.4.6

* Clean up logging. Use --debug to switch it on. Remove deprecated calls to verbose in healpy

## 0.4.3

* Use a common field of view parser with DiSkO
* Add --hdf <filename> option to save the output as an HDF

## 0.4.2

* Read from measurement sets --ms

## 0.4.1

* Use the disko sphere.
* Clean up unused code.
* Use harmonics from disko.
* Specify the --fov and --res as in disko

## 0.4.0

* Move to github repository. Add to pypi. Use disko for utility functions.
* Add a --version CLI argument

## 0.3.0

* Update to python3

## 0.3.3

* Add a gridless binary to plot a GRIDLESS imaging, add a PDF export option

## 0.3.4

* Fix bitrot in multispotless.
