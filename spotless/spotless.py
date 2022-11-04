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
# Tim Molteno 2017-2022 tim@elec.ac.nz
#

import logging
import time
import disko


from scipy.optimize import minimize
import numpy as np
import healpy as hp


from tart.imaging import elaz
from tart.util import constants

from .model import Model
from .source import PointSource

from .sphere import HealpixSphere

logger = logging.getLogger(__name__)
# Add other handlers if you're using this as a library
logger.addHandler(logging.NullHandler())


def get_source_list(source_json, el_limit, jy_limit):
    src_list = []
    if source_json is not None:
        src_list = elaz.from_json(
            source_json, el_limit=el_limit, jy_limit=jy_limit)
    return src_list


class SpotlessBase(object):
    def __init__(self, cal_vis):

        self.vis_arr, baselines = cal_vis.get_all_visibility()
        self.u_arr, self.v_arr, self.w_arr = cal_vis.get_all_uvw()
        logger.info(f"Max u: {np.max(self.u_arr)}")
        self.residual_vis = np.zeros_like(self.vis_arr) + self.vis_arr
        self.model = Model()
        self.working_nside = 2**6

        self.harmonics = {}  # Temporary store for harmonics

    def get_harmonics(self, sphere):
        ''' Create the harmonics for this arrangement of sphere pixels
        '''
        cache_key = sphere.nside
        if (cache_key in self.harmonics):
            return self.harmonics[cache_key]

        n_arr_minus_1 = sphere.n - 1
        harmonic_list = []
        p2j = disko.jomega(constants.L1_FREQ)
        
        for u, v, w in zip(self.u_arr, self.v_arr, self.w_arr):
            harmonic = np.exp(p2j*(u*sphere.l + v*sphere.m + w*n_arr_minus_1))
            harmonic_list.append(harmonic)
        self.harmonics[cache_key] = harmonic_list
        return harmonic_list

    def image_visibilities(self, vis_arr, nside):
        """Create a healpix image from visibilities

        Args:

            vis_arr (np.array): An array of visibilities (the same length as the uvw array
                                in this object
            nside (int):        The healpix nside parameter.
        """

        assert len(vis_arr) == len(self.u_arr)
        logger.info("Imaging Visabilities nside={}".format(nside))
        t0 = time.time()

        sphere = HealpixSphere(nside)

        pixels = np.zeros_like(sphere.visible_pixels)
        harmonic_list = self.get_harmonics(sphere)
        for h, vis in zip(harmonic_list, vis_arr):
            pixels += vis*h
        pixels = np.abs(pixels)
        t1 = time.time()
        logger.info("Seconds elapsed {}".format(t1 - t0))

        sphere.set_visible_pixels(pixels)

        return sphere

    def plot(self, plt, sphere, src_list, show_model):
        rot = (0, 90, 0)
        plt.figure(figsize=(6, 6))
        hp.orthview(sphere.healpix_map, rot=rot, xsize=1000,
                    cbar=True, half_sky=True, hold=True)
        # hp.mollview(healpix_map, rot=rot, xsize=3000, cbar=False)
        hp.graticule(verbose=False)

        if show_model:
            for src in self.model:
                sphere.plot_dot(src.el, src.az)

        if src_list is not None:
            for s in src_list:
                sphere.plot_x(s.el_r, s.az_r)

    def display(self, plt, src_list, nside, show_model):
        sphere = self.image_visibilities(self.residual_vis, nside)
        self.plot(plt, sphere, src_list, show_model)

    def beam(self, plt, nside):
        sphere = self.image_visibilities(
            np.ones_like(self.residual_vis), nside)
        self.plot(plt, sphere, src_list=None, show_model=False)

    def deconvolute(self):
        for i in range(25):
            mod, power, p0 = self.step()
            if power >= p0:
                break
            logger.info("Step {}: Model {}".format(i, mod))
            logger.info("Residual Power {}, dp {}".format(power, p0-power))
        logger.info("Deconvolution Complete")

    def step(self):
        raise NotImplementedError("step is not implemented in the base class")

    @staticmethod
    def power_from_pixels(in_pixels):
        ''' The energy in an image is the mean pixel^2
        '''
        power = np.sum(in_pixels**2)/len(in_pixels)
        return power

    @staticmethod
    def scale_to_power(image_pixels, desired_power):
        ''' Scale the image_pixels so that the pixel_power is in fact
            equal to the desired_power. This is only used for the
            reconstruction process where we use a synthesized beam from
            each source (as the beam profile changes with angle in the sky

            desired_pwer = np.sum((scaling*inpixels)**2)/N
            N*desired_power = np.sum(scaling**2 * inpixels**2)
            N*desired_power = scaling**2 * np.sum(inpixels**2)
            scaling = np.sqrt(N*desired_power / np.sum(inpixels**2)
                    = np.sqrt(desired_power / current_power)
        '''
        current_power = SpotlessBase.power_from_pixels(image_pixels)
        scaling = np.sqrt(desired_power / current_power)
        return image_pixels*scaling

    def pixel_power(self, vis):
        sphere = self.image_visibilities(vis, nside=self.working_nside)
        return sphere.power_from_pixels()

    def vis_power(self, vis):
        '''Calculate the power in an image from a set of visibilities
           This works because the integral of the fourier transform of the
           visibilities can be calculated directly from the visibilities
           without the FT.
           https://courses.cit.cornell.edu/ece531/Lectures/chapter6.pdf

           This is Parseval's Theorem, for the DFT it becomes

           sum_{n=0}^{N-1} abs(x_n)**2 = (1/N) sum_{k=0}{N-1} abs(X_k)**2

           Thus as the visibilities are the F.T of the sky brightness.
           the mean of the abs of the visibilities is proportional
           to the power in the 'image'
        '''
        return np.real(np.sum(vis*np.conj(vis)))

    def power(self, vis):
        return self.vis_power(vis)

    def estimate_initial_point_source(self, vis):
        ''' Estimate an initial point source for the model. This is effectively a brute
            force search over the surface of the sphere
        '''
        sphere = self.image_visibilities(vis, nside=self.working_nside)
        a_0, el_0, az_0 = sphere.get_peak()
        p0 = self.power(vis)
        logger.info(f"Peak: {a_0} vis_pwr={p0} {el_0}")
        return a_0, el_0, az_0, p0

    def get_src_vis(self, src):
        ''' Return the visibility from a single point source
            TODO Incorporate the telescope beam (fringe washing)
        '''
        return src.get_vis(self.u_arr, self.v_arr, self.w_arr)

    def reconstruct_direct(self, nside):
        ''' Reconstruct the image

            Use a thresholded beam from each source as the reconstruction beam. This
            works with non-planar antenna arrays (in which case the beam from each source
            will depend on direction
        '''
        logger.info("Reconstructing Image")

        brightest_source = self.model.brightest()
        weakest_source = self.model.faintest()
        logger.info("Brightest Source {}".format(brightest_source))
        logger.info("Weakest Source {}".format(weakest_source))

        sphere = HealpixSphere(nside)
        harmonic_list = self.get_harmonics(sphere)

        total_source_power = 0.0

        for src in self.model:
            src_vis = self.get_src_vis(src)
            src_map = np.zeros_like(sphere.visible_pixels)

            logger.info("Reconstructing Source {}".format(src))
            for h, vis in zip(harmonic_list, src_vis):
                src_map += vis*h

            # Threshold the src_map to be the central peak
            smap = np.abs(src_map)
            cutoff = np.max(smap)/2
            low_value_flags = smap < cutoff
            smap[low_value_flags] = 0
            sphere.add_visible_pixels(
                SpotlessBase.scale_to_power(smap, src.power))

            total_source_power += src.power

        logger.info("Total Source power {}".format(total_source_power))
        # Smooth the map

        # Get the power
        model_map_power = sphere.power_from_pixels()
        logger.info("model_pixel_map_power {}".format(model_map_power))

        # Add the residual
        residual = self.image_visibilities(self.residual_vis, nside)

        reconstructed_pixels = np.zeros_like(sphere.visible_pixels)
        reconstructed_pixels += residual.visible_pixels

        sphere.add_visible_pixels(np.abs(reconstructed_pixels))

        residual_power = self.power(self.residual_vis)
        logger.info("Residual Power {}".format(residual_power))

        rms_residual = residual.rms()
        peak_reconstructed = np.max(sphere.healpix_map)
        logger.info("RMS residual image {}".format(rms_residual))
        logger.info("Peak reconstructed image {}".format(peak_reconstructed))
        if rms_residual > 0:
            logger.info("SNR {}".format(peak_reconstructed/rms_residual))

        return sphere, model_map_power, residual_power


class Spotless(SpotlessBase):

    def __init__(self, cal_vis):
        super(Spotless, self).__init__(cal_vis)

    def step(self):
        '''
            Return a list of residual visibilities and a source location
            so that the sum of the residual and the source visibilities
            is conserved.
        '''
        d_el = np.radians(1.5)
        a_0, el_0, az_0, p0 = self.estimate_initial_point_source(
            self.residual_vis)
        x0 = [0.1, el_0, az_0]
        logger.info(f"x0 {x0}")

        bounds = []
        src = PointSource(0.1, el_0, az_0)
        for b in src.get_bounds(d_el):
            bounds.append(b)

        fmin = minimize(self.f, x0, method='Nelder-mead', bounds=bounds)
        logger.info(f"fmin {fmin} fun={fmin.fun}")
        a, el, az = fmin.x
        p1 = fmin.fun
        src = PointSource(a, el, az)
        src.power = (p0 - p1)
        logger.info("Source Power {} {}".format(src.power, np.sqrt(a)))
        if src.a > 0.0:
            self.add_source(src)

        return self.model, p1, p0

    def add_source(self, src):
        logger.info("Adding source {}".format(src))
        self.model.add_source(src)
        self.residual_vis -= self.get_src_vis(src)

    def f(self, x):
        '''
            Find the power in the residual
        '''
        src = PointSource(a=x[0], el=x[1], az=x[2])
        pt_vis = self.get_src_vis(src)

        return self.power(self.residual_vis - pt_vis)

    def reconstruct(self, nside):
        logger.info("Reconstructing Image")

        max_u = np.max(self.u_arr)
        max_v = np.max(self.v_arr)
        max_w = np.max(self.w_arr)

        beam_width = np.radians(3.0)  # TODO beamwidth is a function of u,v,w

        brightest_source = self.model.brightest()
        weakest_source = self.model.faintest()
        logger.info("Brightest Source {}".format(brightest_source))
        logger.info("Weakest Source {}".format(weakest_source))

        sphere = HealpixSphere(nside)
        total_source_power = 0.0
        sphere.healpix_map *= 0
        for src in self.model:
            i = sphere.index_of(src.el, src.az)
            sphere.healpix_map[i] += src.power
            total_source_power += src.power

        logger.info("Total Source power {}".format(total_source_power))
        # Smooth the map
        model_pixel_map = hp.sphtfunc.smoothing(
            sphere.healpix_map, fwhm=beam_width, verbose=False)
        # Scale the map so that the total pixel power is correct
        model_pixel_map = np.sqrt(
            np.abs(model_pixel_map))*np.sqrt(len(sphere.visible_pixels))

        # Get the power
        model_map_power = Spotless.power_from_pixels(
            model_pixel_map[sphere.visible_index])
        logger.info("model_pixel_map_power {} {} {}".format(model_map_power,
                                                            total_source_power/model_map_power,
                                                            model_map_power / total_source_power))

        # Add the residual
        residual = self.image_visibilities(self.residual_vis, nside)
        sphere.visible_pixels += residual.visible_pixels

        model_pixel_map[sphere.visible_index] += np.abs(sphere.visible_pixels)

        residual_power = self.power(self.residual_vis)
        logger.info("Residual Power {}".format(residual_power))
        rms_residual = residual.rms()
        peak_reconstructed = np.max(model_pixel_map)
        logger.info("RMS residual image {}".format(rms_residual))
        logger.info("Peak reconstructed image {}".format(peak_reconstructed))
        if rms_residual > 0:
            logger.info("SNR {}".format(peak_reconstructed/rms_residual))

        sphere.healpix_map = model_pixel_map
        return sphere, model_map_power, residual_power
