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


def get_peak(sphere, pixels):
    """The pixel where a positive point source removes the most residual power.

    pixels is a signed back-projection, as SpotlessBase.image_visibilities()
    returns it: its value at a pixel is how well a positive source there matches
    the data, so the largest one is where the next source should start.

    The largest in magnitude will not do. A source that has already been fitted
    leaves negative sidelobes in the residual, and once the remaining positive
    peaks are smaller than those the image's deepest point is a hole: no source
    with a positive amplitude can remove any power there, so the optimizer
    returns nothing and the deconvolution stops with sources still to be found.

    The sphere's own pixels are not signed (image_visibilities() stores their
    magnitude, for display), which is why the image is passed separately.

    Args:
        sphere: The sphere the image was made on.
        pixels: The signed back-projection of the image.

    Returns:
        (value, el, az) of the largest pixel.
    """
    pixels = np.asarray(pixels).ravel()
    i = int(np.argmax(pixels))
    return pixels[i], sphere.el_r[i], sphere.az_r[i]
