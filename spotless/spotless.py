#!/usr/bin/env python
#
# The SPOTLESS algorithm for deconvolution into point sources without gridding.
#
# This relies on a Direct Fourier Transform representation of the image from the visibilities. An
# initial guess is made as to a point source, location amplitude. Then the residual power in the
# image without this source is minimized as a function of the source parameters.
#
# The optimized source is added to the model, and the procedure continues until the residual power
# stops decreasing.
#
# Future Changes:
# Add a gaussial source model (rather than a point source).
#
# Tim Molteno 2017-2026 tim@elec.ac.nz
#

import logging

import healpy as hp
import numpy as np
from disko import DiSkOOperator, vis_to_real
from scipy.optimize import minimize

from .model import Model
from .source import PointSource
from .sphere import get_peak, sphere_copy

logger = logging.getLogger(__name__)


class SpotlessBase(object):
    def __init__(self, disko, sphere, beam=None):
        self.disko = disko
        self.sphere = sphere
        self.el_threshold_r = sphere.el_min_r
        self.vis_arr = np.array(self.disko.vis_arr)
        self.residual_vis = np.zeros_like(self.vis_arr) + self.vis_arr
        self.model = Model()
        self.working_nside = 20
        # An antenna beam pattern (spotless.beam.BeamPattern). None is an
        # antenna that responds equally over the whole sky.
        self.beam = beam

        # The back-projection operator is built on first use and then kept,
        # see _back_projection().
        self._operator = None
        self._operator_arrays = None
        self._operator_frequency = None

    def _back_projection(self, sphere):
        """disko's measurement operator for this sphere, built once and reused.

        Its harmonic blocks depend only on the baselines, the frequency and the
        sphere's pixel geometry, none of which change while a deconvolution
        runs, but DiSkO.image_visibilities() builds a new operator, and so new
        blocks, on every call. This one is rebuilt only when one of those
        inputs is a different object: sphere_copy() shares the geometry arrays,
        so copies of the sphere reuse it. The inputs are kept beside it so that
        an identity cannot be reused by another array.
        """
        arrays = (self.disko.u_arr, self.disko.v_arr, self.disko.w_arr,
                  sphere.l, sphere.m, sphere.n_minus_1, sphere.pixel_areas)
        if (self._operator is None
                or self._operator_frequency != self.disko.frequency
                or any(a is not b for a, b in zip(arrays, self._operator_arrays))):
            self._operator = DiSkOOperator(
                self.disko.u_arr, self.disko.v_arr, self.disko.w_arr,
                self.disko.vis_to_data(), [self.disko.frequency], sphere,
            )
            self._operator_arrays = arrays
            self._operator_frequency = self.disko.frequency
        return self._operator

    def image_visibilities(self, vis_arr, sphere):
        """Create an image from visibilities

        This is the back-projection matched to the model: for the visibilities
        of a PointSource it peaks at that source, which is what
        estimate_initial_point_source() relies on.

        PointSource.get_vis() predicts exp(-j.w.theta). disko's operator is
        the adjoint of a forward model with the conjugate kernel, exp(+j.w.theta),
        so applying it to visibilities v returns sum_i conj(v_i) exp(+j.w.theta_i):
        the conjugate of what is needed, peaking 180 deg in azimuth from a model
        source. Applying it to the conjugated visibilities gives the real part
        of sum_i v_i exp(+j.w.theta_i) instead.

        Args:

            vis_arr (np.array): An array of visibilities (the same length as the uvw array
                                of the disko object
            sphere  (Sphere):   The sphere object that will contain the image.
        """
        assert len(vis_arr) == len(self.disko.u_arr)
        vis = np.conj(vis_arr)
        pixels = self._back_projection(sphere).T @ vis_to_real(vis)
        sphere.set_visible_pixels(np.abs(pixels))
        return pixels.reshape(-1, 1)

    def plot(self, plt, sphere, src_list, show_model):
        sphere.plot(plt, src_list)
        # rot = (0, 90, 0)
        # plt.figure(figsize=(6, 6))
        # hp.orthview(sphere.pixels, rot=rot, xsize=1000,
        #             cbar=True, half_sky=True, hold=True)
        # # hp.mollview(healpix_map, rot=rot, xsize=3000, cbar=False)
        # hp.graticule(verbose=False)
        #
        if show_model:
            for src in self.model:
                sphere.plot_dot(src.el, src.az)

        if src_list is not None:
            for s in src_list:
                sphere.plot_x(plt, s.el_r, s.az_r)

    def display(self, plt, src_list, sphere, show_model):
        self.image_visibilities(self.residual_vis, sphere)
        self.plot(plt, sphere, src_list, show_model)

    def psf(self, plt, sphere):
        """Image the synthesised beam: the response to a unit source at the zenith.

        image_visibilities() replaces the pixels of the sphere it is given, so
        pass a sphere_copy() unless that is wanted.
        """
        self.image_visibilities(np.ones_like(self.residual_vis), sphere)
        self.plot(plt, sphere, src_list=None, show_model=False)

    def beam_pattern(self, plt, sphere):
        """Image the antenna beam pattern over a sphere.

        set_visible_pixels() replaces the pixels of the sphere it is given, so
        pass a sphere_copy() unless that is wanted.
        """
        sphere.set_visible_pixels(self.beam.gain(sphere.el_r, sphere.az_r))
        self.plot(plt, sphere, src_list=None, show_model=False)

    def deconvolute(self, logfile=None, max_steps=50, threshold=0.01):
        """Run the deconvolution.

        Args:
            logfile: Optional file handle to write progress and statistics.
            max_steps: Maximum number of deconvolution steps (default 50).
            threshold: Stop once a step removes less than this fraction of the
                residual power at the start of the run (default 0.01). The source
                fitted in that step stays in the model. 0 leaves stopping to the
                steps that do not improve the residual, and to max_steps.
        """

        def _emit(msg):
            print(msg)
            if logfile is not None:
                logfile.write(msg + "\n")
                logfile.flush()

        stall_count = 0
        patience = 3  # extra steps after first non-improving step
        start_power = self.power(self.residual_vis)  # initial residual power
        best_power = start_power

        for i in range(max_steps):
            mod, power, p0 = self.step()
            if power < best_power:
                # Improvement relative to the best residual power seen.
                stall_count = 0
                dp = best_power - power
                best_power = power
                _emit(
                    f"  Step {i:2d}: residual_power={power:.3f}  removed={dp:.3f}"
                    f"  sources={len(mod.objects)}"
                    f"  nfev={self._opt_nfev}  nit={self._opt_nit}"
                    f"  converged={self._opt_success}"
                )
                if dp < threshold * start_power:
                    _emit(
                        f"  Converged: step {i} removed"
                        f" {100 * dp / start_power:.2f}% of the residual power at the"
                        f" start ({start_power:.3f}), under the"
                        f" {100 * threshold:g}% threshold."
                    )
                    break
            else:
                stall_count += 1
                _emit(
                    f"  Step {i:2d}: residual_power={power:.3f}  no improvement"
                    f"  (stall {stall_count}/{patience})"
                )
                if stall_count >= patience:
                    _emit(
                        f"  Converged: {patience} consecutive steps without"
                        f"  improving on best residual_power={best_power:.3f}."
                    )
                    break
            logger.info("Step {}: Model {}".format(i, mod))
            logger.info("Residual Power {}, dp {}".format(power, p0 - power))
        else:
            _emit(f"  Reached max steps ({max_steps}).")

        for src in self.model:
            src.power = self.vis_power(self.get_src_vis(src))

        # Summarize the final model
        n = len(self.model.objects)
        _emit("")
        if n > 0:
            total = sum(src.power for src in self.model.objects)
            brightest = self.model.brightest()
            faintest = self.model.faintest()
            _emit(f"  Deconvolution complete: {n} sources in model.")
            _emit(f"  Total model power: {total:.2f}")
            _emit(f"  Brightest: a={brightest.a:.3f}  power={brightest.power:.3f}")
            _emit(f"  Faintest:  a={faintest.a:.3f}  power={faintest.power:.3f}")
        else:
            _emit("  Deconvolution complete: no sources found.")
        logger.info("Deconvolution Complete")

    def step(self):
        raise NotImplementedError("step is not implemented in the base class")

    @staticmethod
    def scale_to_power(sphere, desired_power):
        """Scale the sphere_pixels so that the pixel_power is in fact
        equal to the desired_power. This is only used for the
        reconstruction process where we use a synthesized beam from
        each source (as the beam profile changes with angle in the sky

        desired_pwer = np.sum((scaling*inpixels)**2)/N
        N*desired_power = np.sum(scaling**2 * inpixels**2)
        N*desired_power = scaling**2 * np.sum(inpixels**2)
        scaling = np.sqrt(N*desired_power / np.sum(inpixels**2)
                = np.sqrt(desired_power / current_power)
        """
        logger.info(f"scale_to_power({sphere.pixels.shape}, {desired_power})")
        current_power = sphere.get_power()
        if (not np.isfinite(current_power)) or current_power <= 0:
            # An all-zero (or invalid) image has no scale: any factor
            # leaves it unchanged, and dividing by zero would give NaN.
            logger.warning(f"scale_to_power: current power is {current_power}, skipping scaling")
            return
        scaling = np.sqrt(desired_power / current_power)
        sphere.pixels *= scaling

    def pixel_power(self, vis):
        sphere = sphere_copy(self.sphere)
        self.image_visibilities(vis, sphere)
        return sphere.get_power()

    def vis_power(self, vis):
        """Calculate the power in an image from a set of visibilities
        This works because the integral of the fourier transform of the
        visibilities can be calculated directly from the visibilities
        without the FT.
        https://courses.cit.cornell.edu/ece531/Lectures/chapter6.pdf

        This is Parseval's Theorem, for the DFT it becomes

        sum_{n=0}^{N-1} abs(x_n)**2 = (1/N) sum_{k=0}{N-1} abs(X_k)**2

        Thus as the visibilities are the F.T of the sky brightness.
        the sum of the abs of the visibilities is proportional
        to the power in the 'image'
        """
        return np.real(np.sum(vis * np.conj(vis)))

    def power(self, vis):
        return self.vis_power(vis)

    def estimate_initial_point_source(self, vis):
        """Estimate an initial point source for the model. This is effectively a brute
        force search over the surface of the sphere
        """
        sphere = sphere_copy(self.sphere)
        # The signed image: sphere.pixels holds only its magnitude, and the
        # next source has to start at the largest positive pixel, not the
        # largest in magnitude (see get_peak()).
        pixels = self.image_visibilities(vis, sphere)
        a_0, el_0, az_0 = get_peak(sphere, pixels)
        p0 = self.power(vis)
        logger.info(f"Peak: {a_0} vis_pwr={p0} {el_0}")
        return a_0, el_0, az_0, p0

    def get_src_vis(self, src):
        """Return the visibility from a single point source

        Through the antenna beam, if there is one: what a source of flux a*a in
        that direction looks like to the antennas. Every prediction of model
        visibilities comes through here or through model_vis(), so this is the
        one place the beam enters the model. Fringe washing is not modelled.
        """
        vis = src.get_vis(self.disko.u_arr, self.disko.v_arr, self.disko.w_arr)
        if self.beam is not None:
            vis = vis * self.beam.gain(src.el, src.az)
        return vis

    def model_vis(self, model):
        """Return the visibilities of a model's sources, through the beam if
        there is one.

        The same as Model.model_vis(), apart from the beam. MultiSpotless
        refits all the sources together and has to see the beam on every one.
        """
        vis = np.zeros_like(self.disko.u_arr, dtype=np.complex64)
        objects = list(model)
        if self.beam is not None and objects:
            gains = self.beam.gain([src.el for src in objects],
                                   [src.az for src in objects])
        else:
            gains = np.ones(len(objects))
        for src, gain in zip(objects, gains):
            vis += gain * src.get_vis(self.disko.u_arr, self.disko.v_arr,
                                      self.disko.w_arr)
        return vis

    def reconstruct_direct(self):
        """Reconstruct the image

        Use a thresholded beam from each source as the reconstruction beam. This
        works with non-planar antenna arrays (in which case the beam from each source
        will depend on direction
        """
        logger.info("Reconstructing Image")

        brightest_source = self.model.brightest()
        weakest_source = self.model.faintest()
        logger.info("Brightest Source {}".format(brightest_source))
        logger.info("Weakest Source {}".format(weakest_source))

        sphere = sphere_copy(self.sphere)
        total_source_power = 0.0

        scratch = sphere_copy(self.sphere)

        for src in self.model:
            src_vis = self.get_src_vis(src)

            self.image_visibilities(src_vis, scratch)

            # Threshold the scratch to be the central peak
            smap = np.abs(scratch.pixels)
            cutoff = np.max(smap) / 2
            low_value_flags = smap < cutoff
            smap[low_value_flags] = 0

            scratch.pixels = smap
            SpotlessBase.scale_to_power(scratch, src.get_power())
            sphere.pixels += scratch.pixels

            total_source_power += src.get_power()

        del scratch

        logger.info("Total Source power {}".format(total_source_power))

        # Get the power
        model_map_power = sphere.get_power()

        # model_map_power = power_from_pixels(sphere)
        logger.info("model_pixel_map_power {}".format(model_map_power))

        # Add the residual
        residual = sphere_copy(self.sphere)
        self.image_visibilities(self.residual_vis, residual)

        combined = sphere_copy(self.sphere)
        combined.set_visible_pixels(sphere.pixels + residual.pixels, scale=True)

        return sphere, model_map_power, residual.get_power()


class Spotless(SpotlessBase):
    def __init__(self, disko, sphere, beam=None):
        super(Spotless, self).__init__(disko, sphere, beam=beam)

    def step(self):
        """
        Return a list of residual visibilities and a source location
        so that the sum of the residual and the source visibilities
        is conserved.
        """
        d_el = self.disko.get_beam_width().radians() / 2

        a_0, el_0, az_0, p0 = self.estimate_initial_point_source(self.residual_vis)
        a_init = max(a_0, 0.01)
        x0 = [a_init, el_0, az_0]
        logger.info(f"x0 {x0}")

        bounds = []
        src = PointSource(a_init, el_0, az_0)
        for b in src.get_bounds(d_el, self.el_threshold_r):
            logger.info(f"   Bound: {b}")
            bounds.append(b)

        fmin = minimize(
            self.f,
            x0,
            method="L-BFGS-B",
            bounds=bounds,
            options={"ftol": 1e-6, "gtol": 1e-5, "maxiter": 100},
        )
        logger.info(f"fmin {fmin} fun={fmin.fun}")
        self._opt_nfev = fmin.nfev
        self._opt_nit = fmin.nit
        self._opt_success = fmin.success
        a, el, az = fmin.x
        p1 = fmin.fun
        src = PointSource(a, el, az)
        src.power = p0 - p1
        logger.info("Source Power {} {}".format(src.get_power(), np.sqrt(a)))
        if src.a > 0.0:
            self.add_source(src)

        return self.model, p1, p0

    def add_source(self, src):
        logger.info("Adding source {}".format(src))
        self.model.add_source(src)
        self.residual_vis -= self.get_src_vis(src)

    def f(self, x):
        """
        Find the power in the residual
        """
        src = PointSource(a=x[0], el=x[1], az=x[2])
        pt_vis = self.get_src_vis(src)

        return self.power(self.residual_vis - pt_vis)

    def reconstruct(self):
        logger.info("Reconstructing Image")
        sphere = sphere_copy(self.sphere)

        beam_width = self.disko.get_beam_width().radians()
        logger.info(f"Resolution : {beam_width}")

        brightest_source = self.model.brightest()
        weakest_source = self.model.faintest()
        logger.info("Brightest Source {}".format(brightest_source))
        logger.info("Weakest Source {}".format(weakest_source))

        total_source_power = 0.0
        sphere.pixels *= 0
        for src in self.model:
            i = sphere.index_of(src.el, src.az)
            if i < sphere.npix:
                sphere.pixels[i] += src.get_power()
                total_source_power += src.get_power()
            else:
                raise RuntimeError(f"Index {i} >= {sphere.npix}, src={src}")

        logger.info("Total Source power {}".format(total_source_power))

        # Smooth the map
        all_npix = hp.nside2npix(self.sphere.nside)
        all_pixels = np.zeros(all_npix) + hp.UNSEEN
        all_pixels[sphere.pixel_indices] = sphere.pixels

        model_pixel_map = hp.sphtfunc.smoothing(all_pixels, fwhm=beam_width)

        # Scale the map so that the total pixel power is correct
        sphere.pixels = model_pixel_map[sphere.pixel_indices]
        SpotlessBase.scale_to_power(sphere, total_source_power)

        # Get the power
        model_map_power = sphere.get_power()
        logger.info(f"model_pixel_map_power: {model_map_power}")
        logger.info(f"   total_source_power: {total_source_power}")

        # Add the residual
        residual = sphere_copy(self.sphere)
        self.image_visibilities(self.residual_vis, residual)
        residual_power = self.power(self.residual_vis)
        SpotlessBase.scale_to_power(residual, residual_power)
        logger.info(f"       residual_power: {residual_power}")

        combined = sphere_copy(self.sphere)
        combined.pixels = sphere.pixels + residual.pixels
        return combined, model_map_power, residual_power

    def reconstruct_err(self, nside):
        logger.info("Reconstructing Image")

        beam_width = self.disko.get_beam_width().radians()

        brightest_source = self.model.brightest()
        weakest_source = self.model.faintest()
        logger.info("Brightest Source {}".format(brightest_source))
        logger.info("Weakest Source {}".format(weakest_source))

        sphere = sphere_copy(self.sphere)
        total_source_power = 0.0
        sphere.pixels *= 0
        for src in self.model:
            i = sphere.index_of(src.el, src.az)
            sphere.pixels[i] += src.get_power()
            total_source_power += src.get_power()

        logger.info("Total Source power {}".format(total_source_power))

        model_pixel_map = hp.sphtfunc.smoothing(sphere.pixels, fwhm=beam_width)
        # Scale the map so that the total pixel power is correct
        model_pixel_map = np.sqrt(np.abs(model_pixel_map)) * np.sqrt(len(sphere.pixels))

        # Get the power
        model_map_power = Spotless.power_from_pixels(
            model_pixel_map[sphere.pixel_indices]
        )
        logger.info(
            "model_pixel_map_power {} {} {}".format(
                model_map_power,
                total_source_power / model_map_power,
                model_map_power / total_source_power,
            )
        )

        # Add the residual
        residual = sphere_copy(sphere)
        self.image_visibilities(self.residual_vis, residual)
        sphere.pixels += residual.pixels

        # model_pixel_map[sphere.pixel_indices] += np.abs(sphere.pixels)

        residual_power = self.power(self.residual_vis)
        logger.info("Residual Power {}".format(residual_power))
        # rms_residual = rms(residual)
        # peak_reconstructed = np.max(model_pixel_map)
        # logger.info("RMS residual image {}".format(rms_residual))
        # logger.info("Peak reconstructed image {}".format(peak_reconstructed))
        # if rms_residual > 0:
        #     logger.info("SNR {}".format(peak_reconstructed/rms_residual))

        # sphere.pixels = model_pixel_map
        return sphere, model_map_power, residual_power
