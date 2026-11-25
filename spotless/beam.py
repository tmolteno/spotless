#
# Copyright Tim Molteno 2017-2026 tim@elec.ac.nz
# License GPLv3
#
"""Antenna beam patterns for the source model.

A beam pattern is the factor by which an antenna's response scales a source's
flux, as a function of elevation and azimuth. Without one, spotless assumes
every antenna responds equally over the whole sky.

The patterns come from the tart_beam package (github.com/tart-telescope/beams).
--beam takes the same values as tartball's: "tart" for the built-in TART beam,
or the path of a JSON file holding either a beam saved by
tart_beam.Beam.to_json(), or a list of {el, az, gain} samples in degrees, which
is fitted when it is loaded.
"""

import json
import logging

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from tart_beam import Beam, base_tart_beam, elaz_to_vec, fit_from_json

logger = logging.getLogger(__name__)

BUILTIN = ("tart", "default")
"""What --beam takes for the built-in TART beam."""

MIN_GAIN = 0.01
"""The gain is not taken below this above the horizon.

A beam goes to zero at the horizon, and the flux a source is fitted with is
what the data show divided by the gain. So without a floor a noise spike at
the horizon would become a source of enormous flux. Below the horizon the gain
is zero.
"""

GRID_STEP_DEG = 1.0
"""Spacing of the elevation-azimuth table a pattern is looked up in.

Cubic interpolation at this spacing reproduces the built-in TART beam to
within 5e-4 of its response near the horizon, where it curves most, and 6e-6
above 12 degrees elevation. A pattern with structure finer than a few degrees
needs a smaller step.
"""


def load_beam(name):
    """Load the beam pattern that --beam names.

    Args:
        name: "tart" or "default" for the built-in TART beam, or the path of a
            JSON file: a beam saved by tart_beam.Beam.to_json(), or a list of
            {el, az, gain} samples in degrees, which is fitted.

    Returns:
        A BeamPattern.

    Raises:
        ValueError: For anything else, naming it.
    """
    if name in BUILTIN:
        logger.info("Using the built-in TART beam")
        return BeamPattern(base_tart_beam(degree=8, q=1, name="TART_base"),
                           source=name)

    logger.info(f"Loading beam from {name}")
    try:
        with open(name) as f:
            payload = json.load(f)
    except OSError as e:
        raise ValueError(
            f"{name} is neither '{BUILTIN[0]}' nor a readable file ({e.strerror})"
        ) from e
    except ValueError as e:  # json.JSONDecodeError
        raise ValueError(f"{name} is not a beam file: {e}") from e
    try:
        if isinstance(payload, dict) and "coeffs" in payload:
            beam = Beam.from_json(name)
        else:
            beam = fit_from_json(name, degree=8, q=2, ridge=0.0, name="user_beam")
    except (ValueError, KeyError, TypeError, np.linalg.LinAlgError) as e:
        raise ValueError(f"{name} is not a beam file: {e}") from e
    return BeamPattern(beam, source=str(name))


class BeamPattern:
    """A tart_beam.Beam looked up by elevation and azimuth.

    Beam.evaluate() works out its spherical-harmonic basis on every call, which
    costs most of a millisecond for a single direction. The optimiser asks for
    the gain at one new direction on every function evaluation, and
    MultiSpotless does so for every source in its model, tens of thousands of
    times a step. So the beam is evaluated once, over a table of directions,
    and interpolated from that.

    Args:
        beam: The tart_beam.Beam.
        source: What it was loaded from, for the record.
        step_deg: The table's spacing in elevation and in azimuth.
    """

    def __init__(self, beam, source=None, step_deg=GRID_STEP_DEG):
        self.beam = beam
        self.source = source
        el = np.arange(0.0, 90.0 + step_deg / 2, step_deg)
        # The last column repeats the first, so that the table wraps in azimuth.
        az = np.arange(0.0, 360.0 + step_deg / 2, step_deg)
        EL, AZ = np.meshgrid(el, az, indexing="ij")
        table = beam.evaluate(elaz_to_vec(EL.ravel(), AZ.ravel()))
        self._table = RegularGridInterpolator(
            (el, az), table.reshape(EL.shape), method="cubic",
            bounds_error=False, fill_value=0.0,
        )

    def gain(self, el_r, az_r):
        """The beam's response at an elevation and azimuth, in radians.

        Scalars or arrays, broadcast against each other. See MIN_GAIN for what
        happens near and below the horizon.

        The direction goes through tart_beam's own elaz_to_vec(), with the
        elevation and azimuth as spotless measures them: the same call tartball
        makes, so the two read a beam file alike. tart_beam and disko take
        azimuth from different axes, which only shows in a beam that is not
        symmetric about the zenith.
        """
        el_r, az_r = np.broadcast_arrays(np.asarray(el_r, dtype=float),
                                         np.asarray(az_r, dtype=float))
        # disko's azimuths run negative; the table covers 0 to 360 degrees.
        points = np.stack([np.degrees(el_r), np.degrees(az_r) % 360.0], axis=-1)
        single = points.ndim == 1
        gain = self._table(np.expand_dims(points, 0) if single else points)
        gain = np.where(el_r >= 0, np.maximum(gain, MIN_GAIN), 0.0)
        return gain[0] if single else gain

    def to_dict(self):
        """The beam as tart_beam serialises it, and where it came from."""
        ret = self.beam.to_dict()
        ret["source"] = self.source
        return ret

    def __repr__(self):
        return f"BeamPattern({self.beam!r}, source={self.source!r})"
