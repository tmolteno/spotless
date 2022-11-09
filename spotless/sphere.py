#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#
# Classes to hold pixelated spheres
#

import logging
import numpy as np
import healpy as hp

from disko import HealpixSphere

logger = logging.getLogger(__name__)
# Add other handlers if you're using this as a library
logger.addHandler(logging.NullHandler())


def elaz2lmn(el_r, az_r):
    l = np.sin(az_r)*np.cos(el_r)
    m = np.cos(az_r)*np.cos(el_r)
    # Often written in this weird way... np.sqrt(1.0 - l**2 - m**2)
    n = np.sin(el_r)
    return l, m, n


# class Sphere(object):
#     ''' A pixelated Sphere '''
# 
#     def __init__(self, el_r, az_r):
#         self.el_r = el_r
#         self.az_r = az_r
# 
#         self.l, self.m, self.n = elaz2lmn(self.el_r, self.az_r)
# 
#     def get_lmn(self):
#         return self.l, self.m, self.n
# 
#     def index_of(self, el, az):
#         raise NotImplementedError("Use a subclass")
# 
#     @staticmethod
#     def power_from_pixels(_pixels):
#         return np.sum(_pixels**2)/len(_pixels)


# class SpotlessSphere(HealpixSubSphere):
#     ''' A healpix Sphere '''
# 
#     @classmethod
#     def from_resolution(
#         cls, res_arcmin=None, nside=None, theta=0.0, phi=0.0, radius_rad=0.0
#     ):
#         super(SpotlessSphere, cls).from_resolution(res_arcmin, nside, theta, phi, radius_rad)
# 
#     @staticmethod
#     def hp2elaz(theta, phi):
#         el = np.pi/2 - theta
#         az = -phi
#         return el, az
# 
#     @staticmethod
#     def elaz2hp(el, az):
#         theta = np.pi/2 - el
#         phi = -az
#         return theta, phi
# 
#     def plot_dot(self, el, az):
#         theta, phi = self.elaz2hp(el, az)
#         hp.projplot(theta, phi, 'k.', rot=(0, 90, 0))  #
# 
#     def plot_x(self, el, az):
#         theta, phi = self.elaz2hp(el, az)
#         hp.projplot(theta, phi, 'rx', rot=(0, 90, 180))  #
# 
#     def add_pixels(self, pix):
#         self.pixels += pix
# 
def rms(sphere):
    return np.sqrt(np.mean(sphere.pixels**2))

def get_peak(sphere):
    i = np.argmax(sphere.pixels)
    a_0 = sphere.pixels[i]
    el_0 = sphere.el_r[i]
    az_0 = sphere.az_r[i]

    return a_0, el_0, az_0
# 
def power_from_pixels(sphere):
    power = np.sum(sphere.pixels**2)/len(sphere.pixels)
    return power
