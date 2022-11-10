#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#
# Classes to hold pixelated spheres
#

import logging
import numpy as np

logger = logging.getLogger(__name__)
# Add other handlers if you're using this as a library
logger.addHandler(logging.NullHandler())


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
