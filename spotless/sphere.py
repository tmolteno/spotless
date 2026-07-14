#
# Copyright Tim Molteno 2017-2026 tim@elec.ac.nz
# License GPLv3
#
# Classes to hold pixelated spheres
#

import copy
import logging
import numpy as np

logger = logging.getLogger(__name__)


def sphere_copy(sphere):
    """Memory-efficient copy that shares immutable geometry arrays.

    Only the mutable pixels array is duplicated. Coordinate arrays
    (l, m, n, el_r, az_r, pixel_areas, pixel_indices, n_minus_1) are
    shared by reference, reducing per-copy cost from ~9 x npix x 8 bytes
    to npix x 8 bytes.
    """
    ret = copy.copy(sphere)
    ret.pixels = np.array(sphere.pixels)
    return ret


def elaz2lmn(el_r, az_r):
    ll = np.sin(az_r)*np.cos(el_r)
    m = np.cos(az_r)*np.cos(el_r)
    # Often written in this weird way... np.sqrt(1.0 - l**2 - m**2)
    n = np.sin(el_r)
    return ll, m, n


def get_peak(sphere):
    i = np.argmax(sphere.pixels)
    a_0 = sphere.pixels[i]
    el_0 = sphere.el_r[i]
    az_0 = sphere.az_r[i]

    return a_0, el_0, az_0
