# Changes

## 0.9.0

* Fix the deconvolution, which had found (almost) nothing since the move to disko
  1.4.4: `DiSkO.image_visibilities()` became the adjoint of disko's forward model,
  the conjugate of the kernel `PointSource.get_vis()` uses, so the dirty image
  peaked 180 deg in azimuth from every source. `image_visibilities()` now gives it
  the conjugated visibilities, and every image spotless makes goes through it,
  including `reconstruct_direct()`
* Start each source at the brightest pixel of the residual image, not the one
  largest in magnitude. `get_peak()` searched the image as it is displayed,
  magnitudes only, so once a fitted source left negative sidelobes deeper than
  the peaks that remained, every step started at a hole. No source with a
  positive amplitude removes any power there: the optimizer returned nothing,
  and the deconvolution stopped with most of the residual still in it. Every
  subset of the test measurement set tried stopped that way, and nothing else
  stopped a run either, so with the hole gone they ran on to `--max-steps` until
  the threshold below was added. `get_peak()` takes the signed image,
  `get_peak(sphere, pixels)`
* Stop when a step removes less than 1% of the residual power at the start of
  the run. Both variants do, they share the loop; the source fitted in that step
  stays in the model, and `--threshold` sets the percentage (0 turns it off). On
  the test measurement set a deconvolution that stalled on a hole after 10
  sources, with 53% of the power left, now stops after 25 with 29%, and at
  20arcmin with 160 visibilities, where no source was found, after 24 with 11%.
  `--multimodel` on 200 rows of the same set stops after 15 steps in under a
  minute, but on the default 1000 it was still removing over 1% a step at step
  25
* Remove `--show-sources` (from `spotless` and `gridless`), with `--catalog` and
  `--elevation`, which only existed to feed it: sources are drawn by `disko_draw`
  on the file `--HDF` writes. Passing the flag now exits with a pointer to that
  instead of argparse's "unrecognized arguments". This also takes out
  `get_source_list()`, which had been crashing `--file` with a TypeError since
  tart renamed one of its parameters
* Remove `--SVG` for the same reason: `disko_draw` renders SVG from the `--HDF`
  file, and passing it exits with the same pointer. It had not been writing
  where it said either: the files went to the working directory, not `--dir`
* Require disko >= 1.5.2. 1.5.0 raised an AttributeError under pylops 2.9, and
  1.5.2 multiplies the harmonic blocks in float32, 8 to 29 times faster, and no
  longer floors the block size at 256 rows, which made a 2.9 GB transient at
  20arcmin resolution
* Conjugate API data and `--vis` files like `--file` snapshots: they are loaded the
  same way and need it. 9 of the 10 brightest sources in the test snapshot land
  within a beam of a catalogued source with it, 2 without
* Keep one back-projection operator per sphere instead of building one for every
  image, so its harmonic blocks are built once per deconvolution rather than once
  per call. The models are the same bit for bit; the test suite takes 47 s for
  28 tests, from 107 s for 14
* Test the deconvolution itself, which had no tests. A model source is found
  where it was put, alone, as one of two, and with a hole in the residual deeper
  than it, and so is a source injected into the measurement set; the snapshot's
  sources match its catalog, and sources are found in the measurement set. Both
  variants stop at the threshold, the image peaks at a model source, the next
  source starts at its largest positive pixel, `reconstruct_direct()` is centred
  on it, and the operator is reused

## 0.8.0

* Fix #1: `--fits` now fails with a clear "not implemented" message (blocked on
  tmolteno/disko#10) instead of an AttributeError traceback
* Upgrade disko to >=1.4.4: blocked operator caps harmonic cache at 500 MB
* Share immutable geometry arrays across sphere copies (88% less memory per copy)
* Reuse scratch sphere in reconstruct_direct instead of per-source allocation
* Update README memory section to reflect bounded cache

## 0.7.5

* Upgrade disko dependency to >=1.4.4 (blocked operator caps harmonic cache at 500 MB)
* Use shared-geometry sphere copies (88% less memory per copy)
* Reuse scratch sphere in reconstruct_direct (no per-source allocation)
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
