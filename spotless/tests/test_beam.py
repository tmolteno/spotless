#
# Copyright Tim Molteno 2017-2026 tim@elec.ac.nz
# License GPLv3
#
"""Antenna beam patterns: loading them, looking them up, and the model through them."""

import glob
import json
import os
import tempfile
import unittest

import numpy as np
from disko import DiSkO
from disko import HealpixSubFoV as HealpixSubSphere
from tart.operation import settings
from tart_beam import Beam, base_tart_beam, elaz_to_vec
from tart_tools import api_imaging

from spotless import Model, Spotless, source
from spotless.beam import MIN_GAIN, BeamPattern, load_beam
from spotless.multi_spotless import MultiSpotless
from spotless.tests.test_fits_export import run_cli_script
from spotless.tests.test_spotless import SMALL_SPHERE, TEST_JSON, angular_distance

# The CLI on the test snapshot. One step is enough to see what it was given.
CLI = ["--file", TEST_JSON, "--healpix", "--fov", "120deg", "--res", "240arcmin",
       "--max-steps", "1"]


def snapshot_disko():
    """The test snapshot's baselines and visibilities, as spotless_cli loads them."""
    with open(TEST_JSON) as f:
        calib_info = json.load(f)
    config = settings.from_api_json(calib_info["info"]["info"], calib_info["ant_pos"])
    cv, _ = api_imaging.vis_calibrated(
        calib_info["data"][-1][0], config, np.asarray(calib_info["gains"]["gain"]),
        np.asarray(calib_info["gains"]["phase_offset"]), [],
    )
    return DiSkO.from_cal_vis(cv)


class TestLoadBeam(unittest.TestCase):
    """--beam takes the same values as tartball's, and loads them the same way."""

    def setUp(self):
        self.el = np.radians([90.0, 45.0, 5.0])
        self.az = np.radians([0.0, 120.0, 250.0])

    def test_the_built_in_names_give_the_tart_beam(self):
        tart, default = load_beam("tart"), load_beam("default")

        np.testing.assert_array_equal(tart.gain(self.el, self.az),
                                      default.gain(self.el, self.az))
        beam = base_tart_beam(degree=8, q=1, name="TART_base")
        # Looked up from a table, so close to the beam, not equal to it.
        np.testing.assert_allclose(
            tart.gain(self.el, self.az),
            beam.evaluate(elaz_to_vec(np.degrees(self.el), np.degrees(self.az))),
            atol=1e-3,
        )
        self.assertEqual(tart.to_dict()["name"], "TART_base")

    def test_the_tart_beam_tapers_to_the_horizon(self):
        # A least-squares fit to a pattern that is 1 above 10 degrees elevation
        # and falls to 0 at the horizon, so only roughly that.
        pattern = load_beam("tart")

        def gain(el_deg):
            return pattern.gain(np.radians(el_deg), 0.0)

        for el_deg in (90.0, 60.0, 20.0):
            self.assertAlmostEqual(gain(el_deg), 1.0, delta=0.05)
        self.assertAlmostEqual(gain(5.0), 0.5, delta=0.1)
        self.assertEqual(gain(-5.0), 0.0)

    def test_a_saved_beam_is_loaded_as_it_was_saved(self):
        beam = base_tart_beam(degree=4, q=1, name="saved")
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "beam.json")
            beam.to_json(path)
            pattern = load_beam(path)

        np.testing.assert_allclose(
            pattern.gain(self.el, self.az),
            beam.evaluate(elaz_to_vec(np.degrees(self.el), np.degrees(self.az))),
            atol=1e-3,
        )
        self.assertEqual(pattern.source, path)
        self.assertEqual(pattern.to_dict()["name"], "saved")

    def test_sample_records_are_fitted(self):
        # Elevation, azimuth and gain in degrees, as tart_beam reads them.
        records = [{"el": el, "az": az, "gain": 1.0 if el >= 30.0 else 0.0}
                   for el in (90.0, 60.0, 30.0, 0.0)
                   for az in (0.0, 90.0, 180.0, 270.0)]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "samples.json")
            with open(path, "w") as f:
                json.dump(records, f)
            pattern = load_beam(path)

        self.assertEqual(pattern.to_dict()["name"], "user_beam")
        high = pattern.gain(np.radians(60.0), 0.0)
        self.assertGreater(high, 0.5)
        self.assertGreater(high, pattern.gain(np.radians(2.0), 0.0))

    def test_anything_else_is_named_in_the_error(self):
        with self.assertRaisesRegex(ValueError, "no_such_beam.json"):
            load_beam("no_such_beam.json")

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "not_a_beam.json")
            with open(path, "w") as f:
                f.write("[1, 2, 3]")
            with self.assertRaisesRegex(ValueError, "not_a_beam.json"):
                load_beam(path)


class TestBeamPattern(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = load_beam("tart")

    def test_the_lookup_follows_the_beam(self):
        # Interpolated from a table, so close to Beam.evaluate(), not equal.
        # Above 2 degrees, where the floor does not apply, and with disko's
        # negative azimuths.
        rng = np.random.default_rng(1)
        el = np.radians(rng.uniform(2.0, 90.0, 2000))
        az = rng.uniform(-2 * np.pi, 0.0, 2000)

        direct = self.pattern.beam.evaluate(elaz_to_vec(np.degrees(el), np.degrees(az)))

        np.testing.assert_allclose(self.pattern.gain(el, az), direct, atol=1e-3)

    def test_the_gain_is_floored_just_above_the_horizon_and_zero_below(self):
        def gain(el_deg):
            return self.pattern.gain(np.radians(el_deg), 0.0)

        self.assertEqual(gain(-0.5), 0.0)
        self.assertEqual(gain(0.0), MIN_GAIN)
        self.assertEqual(gain(0.5), MIN_GAIN)
        self.assertGreater(gain(3.0), MIN_GAIN)

    def test_azimuth_wraps(self):
        el = np.radians(37.0)
        az = np.radians([-200.0, -1.0, 0.0, 250.0])
        expected = self.pattern.gain(el, az)

        np.testing.assert_allclose(self.pattern.gain(el, az - 2 * np.pi), expected,
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(self.pattern.gain(el, az + 2 * np.pi), expected,
                                   rtol=0, atol=1e-12)

    def test_scalars_and_arrays(self):
        self.assertIsInstance(self.pattern.gain(0.5, 1.0), float)
        self.assertEqual(self.pattern.gain(np.radians([10.0, 20.0, 30.0]), 0.0).shape, (3,))
        self.assertEqual(self.pattern.gain(0.5, np.zeros(4)).shape, (4,))

    def test_directions_go_through_tart_beam_s_own_elaz_to_vec(self):
        # tart_beam and disko take azimuth from different axes. Going through
        # elaz_to_vec() is what makes spotless and tartball read a beam file
        # alike. It only shows in a beam that is not symmetric about the zenith:
        # this one is 1 + 0.49 x, with x = cos(el) cos(az) in tart_beam's frame.
        beam = Beam([(0, 0), (1, 1)], [np.sqrt(4 * np.pi), -1.0], q=0)
        pattern = BeamPattern(beam)
        el, az = np.radians(30.0), np.radians([0.0, 90.0, 180.0])

        gains = pattern.gain(el, az)

        np.testing.assert_allclose(gains, beam.evaluate(elaz_to_vec(30.0, [0.0, 90.0, 180.0])),
                                   atol=1e-4)
        self.assertGreater(gains[0], 1.4)
        self.assertAlmostEqual(gains[1], 1.0, places=3)
        self.assertLess(gains[2], 0.6)

    def test_the_record_says_which_beam_it_was(self):
        record = self.pattern.to_dict()

        self.assertEqual(record["source"], "tart")
        self.assertEqual(record["coeffs"], self.pattern.beam.to_dict()["coeffs"])
        json.dumps(record)  # it goes in the model file


class TestModelThroughBeam(unittest.TestCase):
    """Every prediction of model visibilities goes through the beam."""

    @classmethod
    def setUpClass(cls):
        cls.disko = snapshot_disko()
        cls.beam = load_beam("tart")

    def setUp(self):
        self.sphere = HealpixSubSphere(**SMALL_SPHERE)

    def test_a_source_s_visibilities_are_scaled_by_the_beam(self):
        # 45 degrees is above the taper, 6 degrees is in it.
        for el in (45.0, 6.0):
            with self.subTest(el=el):
                src = source.PointSource(0.5, np.radians(el), np.radians(120.0))

                through = Spotless(self.disko, self.sphere, beam=self.beam).get_src_vis(src)
                plain = Spotless(self.disko, self.sphere).get_src_vis(src)

                np.testing.assert_allclose(through / plain, self.beam.gain(src.el, src.az),
                                           rtol=1e-6)

    def test_multimodel_sees_the_beam_on_every_source(self):
        # MultiSpotless refits all the sources together, so it goes through
        # model_vis() rather than get_src_vis().
        model = Model()
        for el, az in ((6.0, 0.0), (45.0, 90.0), (70.0, 250.0)):
            model.add_source(source.PointSource(0.5, np.radians(el), np.radians(az)))
        u, v, w = self.disko.u_arr, self.disko.v_arr, self.disko.w_arr

        through = MultiSpotless(self.disko, self.sphere, beam=self.beam).model_vis(model)

        expected = np.zeros_like(u, dtype=np.complex64)
        for src in model:
            expected += self.beam.gain(src.el, src.az) * src.get_vis(u, v, w)
        np.testing.assert_allclose(through, expected, atol=1e-6)
        # Without a beam it is what it was.
        np.testing.assert_array_equal(
            MultiSpotless(self.disko, self.sphere).model_vis(model),
            model.model_vis(u, v, w),
        )

    def test_a_source_in_the_taper_gets_its_flux_before_the_beam(self):
        # What the antennas see is a source's flux times the beam. Fitted
        # through the beam it keeps its own flux; fitted without, it loses the
        # gain, a third of it at 6 degrees.
        true = source.PointSource(0.5, np.radians(6.0), np.radians(200.0))
        seen = Spotless(self.disko, self.sphere, beam=self.beam).get_src_vis(true)
        gain = self.beam.gain(true.el, true.az)
        pixel = self.sphere.min_res().radians()

        for with_beam, expected in ((True, true.a), (False, true.a * np.sqrt(gain))):
            with self.subTest(with_beam=with_beam):
                spot = Spotless(self.disko, self.sphere,
                                beam=self.beam if with_beam else None)
                spot.residual_vis = seen.copy()  # a step subtracts from it

                spot.deconvolute()

                brightest = spot.model.brightest()
                self.assertAlmostEqual(brightest.a, expected, delta=0.02 * true.a)
                self.assertLess(
                    angular_distance(brightest.el, brightest.az, true.el, true.az), pixel
                )


class TestBeamOption(unittest.TestCase):
    """--beam and --psf on the command line."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = tmp.name

    def run_spotless(self, *extra, check=True):
        proc = run_cli_script("spotless", "spotless.spotless_cli",
                              CLI + ["--dir", self.dir] + list(extra))
        if check:
            # The deconvolution's own progress goes to stderr as well as to --log.
            self.assertNotIn("Traceback", proc.stderr)
        return proc

    def what_the_model_recorded(self):
        with open(glob.glob(os.path.join(self.dir, "*_model_*.json"))[0]) as f:
            return json.load(f)

    def images(self):
        return [os.path.basename(p) for p in glob.glob(os.path.join(self.dir, "*.png"))]

    def test_bare_beam_means_the_builtin_pattern(self):
        # Writing the images costs a few seconds, so only the first run does.
        for option, images in ((["--beam"], ["--PNG", "--title", "run"]), (["--beam", "tart"], [])):
            with self.subTest(option=option):
                proc = self.run_spotless(*option, *images)

                self.assertEqual(proc.returncode, 0, proc.stderr)
                beam = self.what_the_model_recorded()["beam"]
                self.assertEqual((beam["name"], beam["source"]), ("TART_base", "tart"))
                if images:
                    # The beam image shows the pattern.
                    self.assertTrue(any(n.startswith("run_beam_") for n in self.images()),
                                    self.images())

    def test_a_beam_file_is_loaded(self):
        path = os.path.join(self.dir, "beam.json")
        base_tart_beam(degree=4, q=1, name="from_file").to_json(path)

        proc = self.run_spotless("--beam", path)

        self.assertEqual(proc.returncode, 0, proc.stderr)
        beam = self.what_the_model_recorded()["beam"]
        self.assertEqual((beam["name"], beam["source"]), ("from_file", path))

    def test_a_beam_that_cannot_be_loaded_is_a_usage_error(self):
        proc = self.run_spotless("--beam", "no_such_beam.json", check=False)

        self.assertEqual(proc.returncode, 2)
        self.assertIn("no_such_beam.json", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)

    def test_the_psf_is_imaged(self):
        # This used to be --beam, and it crashed.
        proc = self.run_spotless("--psf", "--PNG", "--title", "run")

        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(any(n.startswith("run_psf_") for n in self.images()), self.images())
        # Without a beam the model does not say anything about one.
        self.assertNotIn("beam", self.what_the_model_recorded())


if __name__ == "__main__":
    unittest.main()
