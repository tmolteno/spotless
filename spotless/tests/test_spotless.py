#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#

import unittest
import logging
import json

import numpy as np

from spotless import source
from spotless import Spotless

from tart.operation import settings
from tart_tools import api_imaging
from tart.imaging import elaz

from disko import DiSkO

logger = logging.getLogger(__name__)
# Add a null handler so logs can go somewhere
logger.addHandler(logging.NullHandler())

TEST_JSON='test_data/test_data.json'


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

        self.spot = Spotless(disko)

    def test_pixel_vis_power(self):
        vis_power = self.spot.vis_power(self.spot.residual_vis)
        pixel_power = self.spot.pixel_power(self.spot.residual_vis)

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
            healpix_map, model_power, residual_power = self.spot.reconstruct(
                nside=nside)
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
            healpix_map, model_power, residual_power = self.spot.reconstruct_direct(
                nside=nside)
            self.assertAlmostEqual(residual_power, 0.0, 2)
            self.assertAlmostEqual(model_power, src_power, 0)
