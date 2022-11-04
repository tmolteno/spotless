#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#

import numpy as np
import logging
import healpy as hp

from scipy.optimize import minimize

from .spotless import SpotlessBase
from .model import Model
from .source import PointSource
from .sphere import HealpixSphere

logger = logging.getLogger(__name__)
# Add other handlers if you're using this as a library
logger.addHandler(logging.NullHandler())


class MultiSpotless(SpotlessBase):

    def __init__(self, cal_vis):
        super(MultiSpotless, self).__init__(cal_vis)

    def step(self):
        '''
            Return a list of residual visibilities and a source location
            so that the sum of the residual and the source visibilities
            is conserved.
        '''
        d_el = np.radians(2.0)

        m_vis = self.model.model_vis(self.u_arr, self.v_arr, self.w_arr)

        a_0, el_0, az_0, p0 = self.estimate_initial_point_source(
            self.vis_arr-m_vis)

        model_vect = self.model.to_vector()
        x0 = np.append(model_vect, [0.1, el_0, az_0])

        bounds = []
        for src in self.model:
            for b in src.get_bounds(d_el):
                bounds.append(b)

        src = PointSource(0.1, el_0, az_0)
        for b in src.get_bounds(d_el):
            bounds.append(b)

        fmin = minimize(self.f_n, x0.flatten(),
                        method='L-BFGS-B', bounds=bounds)
        p1 = fmin.fun

        self.model = Model.from_vector(fmin.x)

        self.residual_vis = self.vis_arr - \
            self.model.model_vis(self.u_arr, self.v_arr, self.w_arr)

        return self.model, p1, p0

    def f_n(self, x):
        '''
            Find the power in the residual
        '''
        mod = Model.from_vector(x)
        m_vis = mod.model_vis(self.u_arr, self.v_arr, self.w_arr)
        return self.power(self.vis_arr - m_vis)

    def reconstruct(self, nside):
        logger.info("Reconstructing Image")
        sphere = HealpixSphere(nside)
        # healpix_map, visible_pixels, el_r, az_r, p2 = sphere.get_pixels()

        max_u = np.max(self.u_arr)
        max_v = np.max(self.v_arr)
        max_w = np.max(self.w_arr)

        beam_width = np.radians(3.0)  # TODO beamwidth is a function of u,v,w

        brightest_source = self.model.brightest()
        weakest_source = self.model.faintest()
        logger.info("Brightest Source {}".format(brightest_source))
        logger.info("Weakest Source {}".format(weakest_source))

        total_source_power = 0.0
        sphere.healpix_map *= 0
        for src in self.model:
            i = sphere.index_of(src.el, src.az)
            src.power = src.a
            sphere.healpix_map[i] += src.a
            total_source_power += src.a

        logger.info("Total Source power {}".format(total_source_power))
        # Smooth the map
        model_pixel_map = hp.sphtfunc.smoothing(
            sphere.healpix_map, fwhm=beam_width, verbose=False)
        # Scale the map so that the total pixel power is correct
        model_pixel_map = np.sqrt(
            np.abs(model_pixel_map))*np.sqrt(len(sphere.visible_pixels))

        # Get the power
        model_map_power = SpotlessBase.power_from_pixels(
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
        logger.info("SNR {}".format(peak_reconstructed/rms_residual))

        sphere.healpix_map = model_pixel_map
        return sphere, model_map_power, residual_power
        # return model_pixel_map, model_map_power, residual_power

    '''
    def old_reconstruct(self, nside):
        logger.info("Reconstructing Image")
        sphere = HealpixSphere(nside)
        healpix_map, visible_pixels, el_r, az_r, p2 = sphere.get_pixels()

        t0 = time.time()

        beam_width = np.radians(3.0)  # TODO beamwidth is a function of u,v,w

        brightest_source = self.model.brightest()
        weakest_source = self.model.faintest()
        logger.info("Brightest Source {}".format(brightest_source))
        logger.info("Weakest Source {}".format(weakest_source))

        healpix_map *= 0
        for src in self.model:
            a = src.a/weakest_source.a

            theta, phi = sphere.elaz2hp(src.el, src.az)
            i = hp.ang2pix(nside, theta, phi)
            healpix_map[i] = a

        healpix_map = hp.sphtfunc.smoothing(
            healpix_map, fwhm=beam_width, verbose=False)
        healpix_map /= np.max(healpix_map)

        # Add the residual weighted so that its peak is equal to the weakest source
        residual, _, _, _ = self.image_visibilities(self.residual_vis, nside)
        residual = residual / np.max(residual)
        visible_pixels += residual*weakest_source

        healpix_map[p2] += np.abs(visible_pixels)

        logger.info("Residual Amplitude {}".format(np.max(residual)))

        return healpix_map, model_map_power, residual_power
    '''
