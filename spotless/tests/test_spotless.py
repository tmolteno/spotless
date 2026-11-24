#
# Copyright Tim Molteno 2017-2026 tim@elec.ac.nz
# License GPLv3
#

import io
import json
import logging
import unittest
from unittest import mock

import numpy as np
from disko import DiSkO
from disko import HealpixSubFoV as HealpixSubSphere
from disko import DiSkOOperator, disko_from_ms
from tart.imaging import elaz
from tart.operation import settings
from tart_tools import api_imaging

from spotless import Spotless, source
from spotless.multi_spotless import MultiSpotless
from spotless.sphere import elaz2lmn, get_peak, sphere_copy

logger = logging.getLogger(__name__)
# Add a null handler so logs can go somewhere
# logger.addHandler(logging.NullHandler())

TEST_JSON = "test_data/test_data.json"
TEST_MS = "test_data/test.ms"

# Coarse enough that deconvolving a whole sky takes a fraction of a second. Sources
# must be well inside it: radius is 85 deg about the zenith.
SMALL_SPHERE = dict(res_arcmin=120, theta=0, phi=0, radius_rad=np.radians(85))


def angular_distance(el_0, az_0, el_1, az_1):
    """Angle (radians) between two directions given as elevation and azimuth."""
    l_0, m_0, n_0 = elaz2lmn(el_0, az_0)
    l_1, m_1, n_1 = elaz2lmn(el_1, az_1)
    return np.arccos(np.clip(l_0 * l_1 + m_0 * m_1 + n_0 * n_1, -1.0, 1.0))


def how_it_ended(log):
    """The line a deconvolution wrote to its log saying why it stopped."""
    lines = [line.strip() for line in log.getvalue().splitlines()
             if line.strip().startswith(("Converged", "Reached"))]
    if len(lines) != 1:
        raise AssertionError(f"expected one line saying why it stopped, got {lines}")
    return lines[0]


class TestSpotless(unittest.TestCase):
    def setUp(self):
        # Load data from a JSON file
        logger.info(f"Getting Data from file: {TEST_JSON}")
        with open(TEST_JSON, "r") as json_file:
            calib_info = json.load(json_file)

        info = calib_info["info"]
        ant_pos = calib_info["ant_pos"]
        config = settings.from_api_json(info["info"], ant_pos)

        flag_list = []

        gains_json = calib_info["gains"]
        gains = np.asarray(gains_json["gain"])
        phase_offsets = np.asarray(gains_json["phase_offset"])
        config = settings.from_api_json(info["info"], ant_pos)

        measurements = []
        for d in calib_info["data"]:
            vis_json, source_json = d
            cv, timestamp = api_imaging.vis_calibrated(
                vis_json, config, gains, phase_offsets, flag_list
            )
            src_list = elaz.from_json(source_json, 0.0)

        disko = DiSkO.from_cal_vis(cv)

        sphere = HealpixSubSphere(
            res_arcmin=60, theta=0, phi=0, radius_rad=np.radians(155)
        )

        self.spot = Spotless(disko, sphere=sphere)

    def test_pixel_vis_power(self):
        vis_power = self.spot.vis_power(self.spot.residual_vis)
        pixel_power = self.spot.pixel_power(self.spot.residual_vis)
        print(f"Vis power {vis_power}, pix {pixel_power}, npix={self.spot.sphere.npix}")
        self.assertAlmostEqual(vis_power / pixel_power, 1.0, 1)

    def test_single_vis_power(self):
        # Point soure at zenith
        src = source.PointSource(0.5, np.pi / 2, 0.0)

        # power in dirty beam of point source
        src_vis = self.spot.get_src_vis(src)
        src_power = self.spot.vis_power(src_vis)

        # Power in the data
        original_power = self.spot.vis_power(self.spot.residual_vis)

        print("Original Power: {}".format(original_power))
        print("Source   Power: {}".format(src_power))

        self.spot.add_source(src)
        new_power = self.spot.vis_power(self.spot.residual_vis)
        print("Updated  Power: {}".format(new_power))

        self.assertTrue(new_power > original_power)
        self.assertAlmostEqual((original_power + src_power) / new_power, 1.0, 1)

    def test_single_image_power(self):
        # Point soure at zenith
        src = source.PointSource(0.5, np.pi / 2, 0.0)

        # power in dirty beam of point source
        src_vis = self.spot.get_src_vis(src)
        src_power = self.spot.pixel_power(src_vis)

        # Power in the data
        original_power = self.spot.pixel_power(self.spot.residual_vis)

        print("Original Power: {}".format(original_power))
        print("Source   Power: {}".format(src_power))

        self.spot.add_source(src)
        new_power = self.spot.pixel_power(self.spot.residual_vis)
        print("Updated  Power: {}".format(new_power))

        # Pixel power is not strictly additive due to cross-terms in the
        # image domain, but should be close.
        self.assertTrue(new_power > original_power)
        self.assertAlmostEqual(original_power + src_power, new_power, delta=2.0)

    def test_single_source_reconstruction(self):
        # Point soure at zenith
        src = source.PointSource(0.5, np.pi / 3, 0.0)
        src_vis = self.spot.get_src_vis(src)
        src_power = self.spot.power(src_vis)
        src.power = src_power
        print("Source   Power: {}".format(src_power))

        self.spot.residual_vis = (
            src_vis  # Create a single noise-free point source image
        )
        residual_power = self.spot.power(self.spot.residual_vis)
        print("Residual Power: {}".format(src_power))

        self.spot.add_source(src)
        residual_power = self.spot.power(self.spot.residual_vis)
        self.assertAlmostEqual(residual_power, 0.0, 2)

        # Now reconstruct
        for nside in [64, 128]:
            sphere = HealpixSubSphere(
                nside=nside, theta=0, phi=0, radius_rad=np.radians(170)
            )
            self.spot.sphere = sphere
            healpix_map, model_power, residual_power = self.spot.reconstruct()
            self.assertAlmostEqual(residual_power, 0.0, 2)
            self.assertAlmostEqual(model_power, src_power, 0)

    def test_single_source_reconstruction_direct(self):
        # Point soure at zenith
        src = source.PointSource(0.5, np.pi / 4, 0.0)
        src_vis = self.spot.get_src_vis(src)
        src_power = self.spot.power(src_vis)
        src.power = src_power
        print("Source   Power: {}".format(src_power))

        self.spot.residual_vis = (
            src_vis  # Create a single noise-free point source image
        )
        residual_power = self.spot.power(self.spot.residual_vis)
        print("Residual Power: {}".format(src_power))

        self.spot.add_source(src)
        residual_power = self.spot.power(self.spot.residual_vis)
        self.assertAlmostEqual(residual_power, 0.0, 2)

        # Now reconstruct
        for nside in [64, 128]:
            sphere = HealpixSubSphere(
                nside=nside, theta=0, phi=0, radius_rad=np.radians(170)
            )
            self.spot.sphere = sphere
            healpix_map, model_power, residual_power = self.spot.reconstruct_direct()
            self.assertAlmostEqual(residual_power, 0.0, 2)
            self.assertAlmostEqual(model_power, src_power, 0)

    def test_image_visibilities_peaks_at_a_model_source(self):
        # The back-projection must be the one matched to PointSource.get_vis():
        # the deconvolution starts each source at the peak of this image. The
        # peak is the pixel holding the source, so it can be a pixel width away.
        self.spot.sphere = HealpixSubSphere(**SMALL_SPHERE)
        src = source.PointSource(0.5, np.pi / 4, np.pi / 4)
        pixel = self.spot.sphere.min_res().radians()

        pixels = self.spot.image_visibilities(self.spot.get_src_vis(src), self.spot.sphere)

        peak = int(np.argmax(pixels.ravel()))
        self.assertLess(
            angular_distance(self.spot.sphere.el_r[peak],
                             self.spot.sphere.az_r[peak], src.el, src.az),
            pixel,
        )

    def test_deconvolution_finds_a_model_source(self):
        """A source made with PointSource.get_vis() is found where it was put.

        disko 1.4.4 changed DiSkO.image_visibilities() to the adjoint of its
        forward model, whose kernel is the conjugate of the one get_vis() uses.
        The dirty image then peaked 180 deg in azimuth from the source, outside
        the bounds step() allows the optimizer, and every step found nothing.
        """
        self.spot.sphere = HealpixSubSphere(**SMALL_SPHERE)
        src = source.PointSource(0.5, np.pi / 4, np.pi / 4)
        self.spot.residual_vis = self.spot.get_src_vis(src)
        original_power = self.spot.power(self.spot.residual_vis)
        beam = self.spot.disko.get_beam_width().radians()

        self.spot.deconvolute()

        # The steps that fail to improve the residual still append sources, but
        # with negligible flux.
        found = [s for s in self.spot.model if s.a > 0.01 * src.a]
        self.assertEqual(len(found), 1)
        self.assertLess(
            angular_distance(found[0].el, found[0].az, src.el, src.az), beam / 10
        )
        self.assertAlmostEqual(found[0].a, src.a, delta=0.01 * src.a)
        self.assertLess(self.spot.power(self.spot.residual_vis),
                        1e-5 * original_power)

    def test_deconvolution_finds_two_model_sources(self):
        # Two sources exercise more than one step: the second must be found in
        # the residual left after the first is subtracted.
        self.spot.sphere = HealpixSubSphere(**SMALL_SPHERE)
        true = [source.PointSource(0.5, np.radians(60), 0.0),
                source.PointSource(np.sqrt(0.125), np.radians(45), np.radians(120))]
        vis = np.zeros_like(self.spot.residual_vis)
        for src in true:
            vis += self.spot.get_src_vis(src)
        self.spot.residual_vis = vis
        beam = self.spot.disko.get_beam_width().radians()

        self.spot.deconvolute()

        found = [s for s in self.spot.model if s.a > 0.1 * true[0].a]
        self.assertEqual(len(found), 2)
        for src in true:
            nearest = min(found, key=lambda s: angular_distance(s.el, s.az, src.el, src.az))
            self.assertLess(
                angular_distance(nearest.el, nearest.az, src.el, src.az), beam / 10
            )
            self.assertAlmostEqual(nearest.a, src.a, delta=0.05 * src.a)

    def test_get_peak_is_the_largest_positive_pixel(self):
        # image_visibilities() returns a signed image, as a column, and stores
        # only its magnitude in the sphere.
        sphere = HealpixSubSphere(**SMALL_SPHERE)
        pixels = np.zeros((sphere.npix, 1))
        hole, peak = 10, 20
        pixels[hole] = -100.0
        pixels[peak] = 3.0

        a_0, el_0, az_0 = get_peak(sphere, pixels)

        self.assertEqual(a_0, 3.0)
        self.assertEqual(el_0, sphere.el_r[peak])
        self.assertEqual(az_0, sphere.az_r[peak])

    def test_a_hole_does_not_hide_the_next_source(self):
        """The next source starts at the largest positive pixel of the image.

        image_visibilities() stores the magnitude of the back-projection in the
        sphere, for display, and get_peak() used to search that. A source that
        has been subtracted leaves negative sidelobes in the residual, and once
        the deepest of them is deeper than the peaks that remain, every step
        started at a hole. No source with a positive amplitude removes any power
        there, so the optimizer returned nothing and the deconvolution stopped
        with a source still in the residual. Every subset of the test
        measurement set tried stopped that way, most of the power still in it.
        """
        self.spot.sphere = HealpixSubSphere(**SMALL_SPHERE)
        # A residual with a source in it and a hole deeper than it: a source
        # has been subtracted, leaving its negative beam, and another is still
        # to be found.
        wanted = source.PointSource(0.5, np.radians(60), 0.0)
        subtracted = source.PointSource(0.6, np.radians(45), np.radians(120))
        self.spot.residual_vis = (self.spot.get_src_vis(wanted)
                                  - self.spot.get_src_vis(subtracted))
        original_power = self.spot.power(self.spot.residual_vis)
        beam = self.spot.disko.get_beam_width().radians()
        pixels = self.spot.image_visibilities(self.spot.residual_vis,
                                              sphere_copy(self.spot.sphere))
        # The case only arises when the hole is the deepest point of the image.
        self.assertGreater(-np.min(pixels), np.max(pixels))

        self.spot.deconvolute()

        brightest = self.spot.model.brightest()
        self.assertIsNotNone(brightest, "no source was found")
        self.assertLess(
            angular_distance(brightest.el, brightest.az, wanted.el, wanted.az), beam / 10
        )
        self.assertAlmostEqual(brightest.a, wanted.a, delta=0.05 * wanted.a)
        # The source is a third of the residual power. The steps after it add
        # faint sources around the hole, which a positive amplitude cannot
        # model, and remove some more.
        self.assertLess(self.spot.power(self.spot.residual_vis), 0.8 * original_power)

    def test_stops_when_a_step_removes_under_the_threshold(self):
        """The run ends at the first step to remove less than 1% of the residual
        power it started with.

        Once the sources are fitted, later steps only fit what the earlier fits
        left behind, and on real data that went on to --max-steps.
        """
        self.spot.sphere = HealpixSubSphere(**SMALL_SPHERE)
        # The third source is 0.2% of the starting power: it, or anything else
        # left once the first two are fitted, is under the threshold.
        true = [source.PointSource(0.5, np.radians(60), 0.0),
                source.PointSource(0.2, np.radians(45), np.radians(120)),
                source.PointSource(np.sqrt(0.0005), np.radians(70), np.radians(240))]
        self.spot.residual_vis = sum(self.spot.get_src_vis(src) for src in true)
        beam = self.spot.disko.get_beam_width().radians()

        log = io.StringIO()
        self.spot.deconvolute(logfile=log)

        self.assertIn("under the 1% threshold", how_it_ended(log))
        # The step that stopped it stays in the model, and nothing follows it.
        self.assertEqual(len(self.spot.model.objects), 3)
        for src in true[:2]:
            nearest = min(self.spot.model.objects,
                          key=lambda s: angular_distance(s.el, s.az, src.el, src.az))
            self.assertLess(
                angular_distance(nearest.el, nearest.az, src.el, src.az), beam / 10
            )
            self.assertAlmostEqual(nearest.a, src.a, delta=0.05 * src.a)

        # 0 leaves stopping to the steps that do not improve the residual, and
        # to max_steps.
        fresh = Spotless(self.spot.disko, self.spot.sphere)
        fresh.residual_vis = sum(fresh.get_src_vis(src) for src in true)
        log = io.StringIO()
        fresh.deconvolute(logfile=log, max_steps=6, threshold=0)

        self.assertNotIn("threshold", how_it_ended(log))
        self.assertGreater(len(fresh.model.objects), 3)

    def test_reconstruct_direct_is_centred_on_the_source(self):
        # Each source's thresholded beam is placed where image_visibilities()
        # puts it, so a conjugated kernel scattered the reconstruction 180 deg
        # in azimuth from every source.
        self.spot.sphere = HealpixSubSphere(**SMALL_SPHERE)
        src = source.PointSource(0.5, np.pi / 4, np.pi / 4)
        src_vis = self.spot.get_src_vis(src)
        src.power = self.spot.power(src_vis)
        self.spot.residual_vis = src_vis
        self.spot.add_source(src)
        pixel = self.spot.sphere.min_res().radians()

        image, _, _ = self.spot.reconstruct_direct()

        peak = int(np.argmax(image.pixels))
        self.assertLess(
            angular_distance(image.el_r[peak], image.az_r[peak], src.el, src.az),
            pixel,
        )

    def test_snapshot_sources_match_the_catalog(self):
        """Deconvolving the test snapshot finds catalogued sources where they are.

        This checks the whole convention against the catalog saved with the
        data: the DiSkO.from_cal_vis() baselines, the conjugation spotless_cli
        applies to them, the back-projection and PointSource.get_vis(). Getting
        any of them wrong puts the sources 180 deg off in azimuth.
        """
        self.spot.sphere = HealpixSubSphere(**SMALL_SPHERE)
        # What spotless_cli does to data loaded with DiSkO.from_cal_vis().
        self.spot.residual_vis = np.conj(self.spot.residual_vis)
        with open(TEST_JSON, "r") as json_file:
            catalog = [s for s in json.load(json_file)["data"][-1][1] if s["el"] > 20.0]
        beam = self.spot.disko.get_beam_width().radians()

        self.spot.deconvolute()

        brightest = sorted(self.spot.model.objects, key=lambda s: -s.power)[:10]
        matched = sum(
            any(angular_distance(s.el, s.az, np.radians(c["el"]), np.radians(c["az"])) < beam
                for c in catalog)
            for s in brightest
        )
        self.assertGreaterEqual(matched, 8)

    def test_back_projection_operator_is_reused(self):
        # The harmonic blocks are the expensive part of imaging. They must be
        # built once for a sphere and its copies, not once per image, and
        # rebuilt when the geometry they depend on is replaced.
        self.spot.sphere = HealpixSubSphere(**SMALL_SPHERE)
        with mock.patch("spotless.spotless.DiSkOOperator", wraps=DiSkOOperator) as build:
            for _ in range(3):
                self.spot.image_visibilities(self.spot.residual_vis,
                                             sphere_copy(self.spot.sphere))
            self.assertEqual(build.call_count, 1)

            self.spot.sphere = HealpixSubSphere(**SMALL_SPHERE)
            self.spot.image_visibilities(self.spot.residual_vis, self.spot.sphere)
            self.assertEqual(build.call_count, 2)


class TestMeasurementSet(unittest.TestCase):
    """Deconvolve the test measurement set.

    Unlike a snapshot, a measurement set is loaded as explicit conjugate
    pairs. disko_from_ms() draws its visibilities at random unless it is given
    a generator, and what the deconvolution finds depends a lot on the draw, so
    one test injects a source that dominates the data and the other fixes it.
    """

    def setUp(self):
        self.sphere = HealpixSubSphere(**SMALL_SPHERE)

    def _disko(self, num_vis, rng):
        return disko_from_ms(TEST_MS, column="DATA", num_vis=num_vis,
                             res=self.sphere.min_res(), channel=0, field_id=0,
                             ddid=0, rng=rng)

    def test_an_injected_source_is_found_where_it_was_put(self):
        # Four times the power of the whole measurement set, so the draw
        # doesn't matter. get_vis() scales with a**2, so the power is n * a**4.
        spot = Spotless(self._disko(200, None), self.sphere)
        n = len(spot.residual_vis)
        src = source.PointSource(
            (4 * spot.power(spot.residual_vis) / n) ** 0.25,
            np.radians(45), np.radians(90),
        )
        spot.residual_vis = spot.residual_vis + spot.get_src_vis(src)
        beam = spot.disko.get_beam_width().radians()

        spot.deconvolute()

        brightest = max(spot.model.objects, key=lambda s: s.a)
        self.assertLess(
            angular_distance(brightest.el, brightest.az, src.el, src.az), beam / 5
        )
        self.assertAlmostEqual(brightest.a, src.a, delta=0.1 * src.a)

    def test_sources_are_found(self):
        spot = Spotless(self._disko(200, np.random.default_rng(2)), self.sphere)
        original_power = spot.power(spot.residual_vis)

        spot.deconvolute()

        brightest = max(s.power for s in spot.model.objects)
        found = [s for s in spot.model.objects if s.power > 0.005 * brightest]
        self.assertGreaterEqual(len(found), 5)
        # Starting steps at holes instead of positive peaks used to stop this
        # draw after 10 sources, with over half the power still in the residual.
        self.assertLess(spot.power(spot.residual_vis), 0.4 * original_power)
        for s in found:
            self.assertGreater(s.el, 0.0)

    def test_multimodel_stops_at_the_threshold(self):
        # MultiSpotless stops the same way: it runs the same deconvolute(). Its
        # step is scripted out here, because on a measurement set it takes
        # minutes to get to the threshold, and where it gets there depends on
        # the draw.
        spot = MultiSpotless(self._disko(80, np.random.default_rng(2)), self.sphere)
        start = spot.power(spot.residual_vis)
        left = [0.9, 0.8, 0.795]  # the residual power after each step

        def scripted_step():
            spot.model.add_source(source.PointSource(0.1, 0.5, 0.5))
            spot._opt_nfev = spot._opt_nit = 0  # the log line prints them
            spot._opt_success = True
            return spot.model, left.pop(0) * start, start

        log = io.StringIO()
        with mock.patch.object(spot, "step", scripted_step):
            spot.deconvolute(logfile=log, max_steps=10)

        self.assertIn("under the 1% threshold", how_it_ended(log))
        self.assertEqual(len(spot.model.objects), 3)
