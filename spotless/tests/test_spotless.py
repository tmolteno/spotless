#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#

import json
import logging
import unittest

import numpy as np
from disko import DiSkO, HealpixSubSphere
from tart.imaging import elaz
from tart.operation import settings
from tart_tools import api_imaging

from spotless import Spotless, source

logger = logging.getLogger(__name__)
# Add a null handler so logs can go somewhere
# logger.addHandler(logging.NullHandler())

TEST_JSON = 'test_data/test_data.json'


class TestSpotless(unittest.TestCase):

    def setUp(self):
        # Load data from a JSON file
        logger.info(f"Getting Data from file: {TEST_JSON}")
        with open(TEST_JSON, 'r') as json_file:
            calib_info = json.load(json_file)

        info = calib_info['info']
        ant_pos = calib_info['ant_pos']
        config = settings.from_api_json(info['info'], ant_pos)

        flag_list = []

        gains_json = calib_info['gains']
        gains = np.asarray(gains_json['gain'])
        phase_offsets = np.asarray(gains_json['phase_offset'])
        config = settings.from_api_json(info['info'], ant_pos)

        measurements = []
        for d in calib_info['data']:
            vis_json, source_json = d
            cv, timestamp = api_imaging.vis_calibrated(
                vis_json, config, gains, phase_offsets, flag_list)
            src_list = elaz.from_json(source_json, 0.0)

        disko = DiSkO.from_cal_vis(cv)

        sphere = HealpixSubSphere(
            res_arcmin=60, theta=0, phi=0, radius_rad=np.radians(155))

        self.spot = Spotless(disko, sphere=sphere)

    def test_pixel_vis_power(self):
        vis_power = self.spot.vis_power(self.spot.residual_vis)
        pixel_power = self.spot.pixel_power(self.spot.residual_vis)
        print(f"Vis power {vis_power}, pix {pixel_power}, npix={self.spot.sphere.npix}")
        self.assertAlmostEqual(vis_power/pixel_power, 1.0, 1)

    def test_single_vis_power(self):
        # Point soure at zenith
        src = source.PointSource(0.5, np.pi/2, 0.0)

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
        self.assertAlmostEqual((original_power + src_power)/new_power, 1.0, 1)

    def test_single_image_power(self):
        # Point soure at zenith
        src = source.PointSource(0.5, np.pi/2, 0.0)

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

        test_power = original_power + src_power

        self.assertTrue(new_power > original_power)
        self.assertAlmostEqual(original_power + src_power, new_power, 0)

    def test_single_source_reconstruction(self):
        # Point soure at zenith
        src = source.PointSource(0.5, np.pi/3, 0.0)
        src_vis = self.spot.get_src_vis(src)
        src_power = self.spot.power(src_vis)
        src.power = src_power
        print("Source   Power: {}".format(src_power))

        self.spot.residual_vis = src_vis  # Create a single noise-free point source image
        residual_power = self.spot.power(self.spot.residual_vis)
        print("Residual Power: {}".format(src_power))

        self.spot.add_source(src)
        residual_power = self.spot.power(self.spot.residual_vis)
        self.assertAlmostEqual(residual_power, 0.0, 2)

        # Now reconstruct
        for nside in [64, 128]:
            sphere = HealpixSubSphere(nside=nside,
                                      theta=0, phi=0,
                                      radius_rad=np.radians(170))
            self.spot.sphere = sphere
            healpix_map, model_power, residual_power = self.spot.reconstruct()
            self.assertAlmostEqual(residual_power, 0.0, 2)
            self.assertAlmostEqual(model_power, src_power, 0)

    def test_single_source_reconstruction_direct(self):
        # Point soure at zenith
        src = source.PointSource(0.5, np.pi/4, 0.0)
        src_vis = self.spot.get_src_vis(src)
        src_power = self.spot.power(src_vis)
        src.power = src_power
        print("Source   Power: {}".format(src_power))

        self.spot.residual_vis = src_vis  # Create a single noise-free point source image
        residual_power = self.spot.power(self.spot.residual_vis)
        print("Residual Power: {}".format(src_power))

        self.spot.add_source(src)
        residual_power = self.spot.power(self.spot.residual_vis)
        self.assertAlmostEqual(residual_power, 0.0, 2)

        # Now reconstruct
        for nside in [64, 128]:
            sphere = HealpixSubSphere(nside=nside,
                                      theta=0, phi=0,
                                      radius_rad=np.radians(170))
            self.spot.sphere = sphere
            healpix_map, model_power, residual_power = self.spot.reconstruct_direct()
            self.assertAlmostEqual(residual_power, 0.0, 2)
            self.assertAlmostEqual(model_power, src_power, 0)
