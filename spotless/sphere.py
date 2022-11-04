#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#
# Classes to hold pixelated spheres
#

import logging
import numpy as np
import healpy as hp

logger = logging.getLogger(__name__)
# Add other handlers if you're using this as a library
logger.addHandler(logging.NullHandler())


def elaz2lmn(el_r, az_r):
    l = np.sin(az_r)*np.cos(el_r)
    m = np.cos(az_r)*np.cos(el_r)
    # Often written in this weird way... np.sqrt(1.0 - l**2 - m**2)
    n = np.sin(el_r)
    return l, m, n


class Sphere(object):
    ''' A pixelated Sphere '''

    def __init__(self, el_r, az_r):
        self.el_r = el_r
        self.az_r = az_r

        self.l, self.m, self.n = elaz2lmn(self.el_r, self.az_r)

    def get_lmn(self):
        return self.l, self.m, self.n

    def index_of(self, el, az):
        raise NotImplementedError("Use a subclass")

    @staticmethod
    def power_from_pixels(_pixels):
        return np.sum(_pixels**2)/len(_pixels)


class HealpixSphere(Sphere):
    ''' A healpix Sphere '''

    def __init__(self, nside):
        self.nside = nside
        self.npix = hp.nside2npix(self.nside)

        all_pixels = np.arange(self.npix)
        theta, phi = hp.pix2ang(nside, all_pixels)
        logger.info(f"theta {theta}")

        # Find all the pixels above the horizon
        self.visible_index = np.flatnonzero(theta < np.pi/2)

        self.healpix_map = np.zeros(self.npix)  # + hp.UNSEEN
        self.visible_pixels = np.zeros_like(
            self.visible_index, dtype=np.complex64)
        el_r, az_r = self.hp2elaz(
            theta[self.visible_index], phi[self.visible_index])
        logger.info(f"self.visible_index {self.visible_index}")
        logger.info(f"theta {theta[self.visible_index]}")

        super(HealpixSphere, self).__init__(el_r, az_r)

    @staticmethod
    def hp2elaz(theta, phi):
        el = np.pi/2 - theta
        az = -phi
        return el, az

    @staticmethod
    def elaz2hp(el, az):
        theta = np.pi/2 - el
        phi = -az
        return theta, phi

    def index_of(self, el, az):
        theta, phi = self.elaz2hp(el, az)
        return hp.ang2pix(self.nside, theta, phi)

    def plot_dot(self, el, az):
        theta, phi = self.elaz2hp(el, az)
        hp.projplot(theta, phi, 'k.', rot=(0, 90, 0))  #

    def plot_x(self, el, az):
        theta, phi = self.elaz2hp(el, az)
        hp.projplot(theta, phi, 'rx', rot=(0, 90, 180))  #

    def set_visible_pixels(self, pix):
        self.healpix_map[self.visible_index] = pix

    def add_visible_pixels(self, pix):
        self.healpix_map[self.visible_index] += pix

    def rms(self):
        _pixels = self.healpix_map[self.visible_index]
        return np.sqrt(np.mean(_pixels**2))

    def get_peak(self):
        i = np.argmax(self.healpix_map[self.visible_index])
        a_0 = self.healpix_map[i]
        el_0 = self.el_r[i]
        az_0 = self.az_r[i]

        return a_0, el_0, az_0

    def power_from_pixels(self):
        _pixels = self.healpix_map[self.visible_index]
        # super(HealpixSphere, self).power_from_pixels(_pixels)
        power = np.sum(_pixels**2)/len(_pixels)
        return power
