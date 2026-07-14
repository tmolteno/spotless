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

## Command Line Usage

### Data Sources

Spotless can read visibilities from three sources:

| Source | Flag | Example |
|--------|------|---------|
| TART API | `--api` | `spotless --api https://tart.elec.ac.nz/signal` |
| CASA Measurement Set | `--ms` | `spotless --ms test_data/test.ms` |
| JSON snapshot file | `--file` | `spotless --file observation.json` |

### Imaging Options

```
  --fov FOV             Field of view (e.g., 160deg)
  --res RES             Resolution (e.g., 120arcmin)
  --healpix             Use HEALPix pixelisation
  --nvis NVIS           Number of visibilities to use (default: 1000)
  --channel CHANNEL     Frequency channel (default: 0)
```

### Output Formats

```
  --display             Show image interactively
  --PNG                 Save as PNG
  --SVG                 Save as SVG
  --PDF                 Save as PDF
  --fits                Save as FITS
  --HDF FILENAME        Save field of view as HDF5
  --save-model-json FILE  Save point-source model as JSON
  --dir DIR             Output directory (default: .)
  --title TITLE         Prefix for output filenames
```

### Algorithms

| Flag | Description |
|------|-------------|
| (default) | Sequential Spotless — finds one source at a time |
| `--multimodel` | MultiSpotless — jointly optimises all sources |

### Other Flags

```
  --show-sources        Overlay known sources from catalog
  --show-model          Show the model source locations
  --elevation ELEV      Elevation limit for source display (degrees, default: 20)
  --beam                Generate a dirty beam image
  --max-steps N         Maximum deconvolution steps (default: 50)
  --log FILE            Save deconvolution statistics to FILE
  --version             Print version and exit
  --debug               Enable debug logging
```

### Examples

```sh
# Image a measurement set with MultiSpotless, save as SVG
spotless --ms test_data/test.ms --healpix --fov 160deg --res 120arcmin \
    --multimodel --SVG --title my_source

# Live all-sky imaging from the TART telescope
spotless --api https://tart.elec.ac.nz/signal --display --show-sources \
    --healpix --fov 180deg --res 60arcmin

# Calibrate using the spotless model
spotless_calibrate --api https://tart.elec.ac.nz/signal
```

## Performance

### Compute

Spotless benefits from multi-threaded BLAS libraries for numpy operations.
Set the following environment variables to use multiple cores:

```sh
# Linux / macOS
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4

# Or set for a single run
OPENBLAS_NUM_THREADS=4 spotless --ms data.ms --healpix
```

For large measurement sets, use `--nvis` to limit the number of
visibilities used in the peak search (the optimizer uses all
visibilities for accuracy).

### Memory

Since disko 1.4.4, the harmonic cache uses a blocked, matrix-free operator
with a hard 500 MB cap.  Memory no longer scales as O(n_vis x n_pix);
instead it is bounded by the block cache plus per-sphere pixel arrays.

Sphere copies share immutable geometry arrays (l, m, n, el_r, az_r,
pixel_areas) by reference, so each copy costs only `n_pix x 8` bytes
(the pixels array) instead of ~9 x `n_pix x 8` bytes.

| nside | n_pix | per sphere copy | harmonic cache |
|-------|-------|----------------|----------------|
| 64 | 49,152 | 0.4 MB | <= 500 MB |
| 128 | 196,608 | 1.6 MB | <= 500 MB |
| 256 | 786,432 | 6.3 MB | <= 500 MB |
| 512 | 3,145,728 | 25.2 MB | <= 500 MB |

Guidelines:

- Use the coarsest resolution acceptable for your science (`--res`).
- Limit visibilities with `--nvis` for initial exploration.
- Use `--fov` to restrict the field of view for targeted high-resolution imaging.

## Documentation

A LaTeX article describing the algorithm in detail is available in `doc/spotless.tex`.
Build with:

    cd doc && make

Generate example images and convert to PNG:

    cd doc && make images && make pngs

## TODO

* Add Gaussian Source Model
* Make explicit the antenna model (gain as a function of angular coordinates). We are assuming it is hemispherical here.
* Prove the relationship between power in the image and visibilty amplitudes. This might only work when the image 
  tends towards a random one. But this is OK since as we remove the sources the residual becomes more and more random.
* Run an MCMC on the multimodel option to estimate uncertainty in the model. Then use this uncertainty as a stopping criterion (when new model components no longer have certain amplitude or position)

## Author

* Tim Molteno (tim@elec.ac.nz)

## Development work
    
If you are developing this package, install uv and then:
```
	make sync
```
This installs all dependencies (including dev tools like flake8) in a virtual environment.
To run tests: `make test`
To lint: `make lint`

## Changes

See [CHANGES.md](CHANGES.md).
