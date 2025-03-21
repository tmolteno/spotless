#!/usr/bin/env python
import matplotlib
import os
if os.name == 'posix' and "DISPLAY" not in os.environ:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
import json
import logging
from copy import deepcopy

import numpy as np

from tart.operation import settings

from tart_tools import api_handler
from tart_tools import api_imaging
from tart.imaging import elaz

from .spotless import Spotless
from .spotless import get_source_list


logger = logging.getLogger()


def handle_image(args, img, title, time_repr, src_list=None):
    """ This function manages the output of an image, drawing sources e.t.c."""
    image_title = '{}_{}'.format(title, time_repr)
    plt.title(image_title)
    if args.fits:
        fname = '{}.fits'.format(image_title)
        fpath = os.path.join(args.dir, fname)
        api_imaging.save_fits_image(img, fname=fname, out_dir=args.dir, timestamp=time_repr)
        logger.info("Generating {}".format(fname))
    if args.PNG:
        fname = '{}.png'.format(image_title)
        fpath = os.path.join(args.dir, fname)
        plt.savefig(fpath, dpi=300)
        logger.info("Generating {}".format(fname))
    if args.PDF:
        fname = '{}.pdf'.format(image_title)
        fpath = os.path.join(args.dir, fname)
        plt.savefig(fpath, dpi=600)
        logger.info("Generating {}".format(fname))
    if args.display:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Generate an all-sky TART Gridless Image',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '--api', required=False,
        default='https://tart.elec.ac.nz/signal',
        help="Telescope API server URL.")

    parser.add_argument('--catalog', required=False, default='https://tart.elec.ac.nz/catalog', help="Catalog API URL.")

    parser.add_argument('--file', required=False, default=None, help="Snapshot ovservation saved JSON file (visiblities, positions and more).")
    parser.add_argument('--vis', required=False, default=None, help="Use a local JSON file containing the visibilities to create the image.")
    parser.add_argument('--dir', required=False, default='.', help="Output directory.")
    parser.add_argument('--rotation', type=float, default=0.0, help="Apply rotation (in degrees) to the antenna positions.")
    parser.add_argument('--nside', type=int, default=8, help="Healpix nside parameter for display purposes only.")

    parser.add_argument('--beam', action="store_true", help="Generate a gridless beam.")

    parser.add_argument('--elevation', type=float, default=20.0, help="Elevation limit for displaying sources (degrees)")
    parser.add_argument('--display', action="store_true", help="Display Image to the user")
    parser.add_argument('--fits', action="store_true", help="Generate a FITS format image")
    parser.add_argument('--PNG', action="store_true", help="Generate a PNG format image")
    parser.add_argument('--PDF', action="store_true", help="Generate a PDF format image")
    parser.add_argument('--show-sources', action="store_true", help="Show known sources on images (only works on PNG).")

    parser.add_argument('--title', required=False, default="", help="Prefix the title.")

    parser.add_argument('--version', action="store_true",
                        help="Display the current version")
    parser.add_argument('--debug', action="store_true",
                        help="Display debugging information")

    source_json = None

    ARGS = parser.parse_args()

    if ARGS.debug:
        level = logging.DEBUG
    else:
        level = logging.ERROR

    logger = logging.getLogger('spotless')
    logger.setLevel(level)

    if ARGS.debug:
        fh = logging.FileHandler('spotless.log')
        fh.setLevel(level)

        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(level)

        # create formatter
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        # add formatter to ch
        ch.setFormatter(formatter)
        fh.setFormatter(formatter)

        # add ch to logger
        logger.addHandler(ch)

    if ARGS.version:
        version = version("spotless")
        print(f"gridless: Version {version}")
        print("          (c) 2022-2023 Tim Molteno")
        sys.exit(0)

    
    if ARGS.file:
        logger.info("Getting Data from file: {}".format(ARGS.file))
        # Load data from a JSON file
        with open(ARGS.file, 'r') as json_file:
            calib_info = json.load(json_file)

        info = calib_info['info']
        ant_pos = calib_info['ant_pos']
        config = settings.from_api_json(info['info'], ant_pos)

        flag_list = [] # [4, 5, 14, 22]

        original_positions = deepcopy(config.get_antenna_positions())

        gains_json = calib_info['gains']
        gains = np.asarray(gains_json['gain'])
        phase_offsets = np.asarray(gains_json['phase_offset'])
        config = settings.from_api_json(info['info'], ant_pos)
    
        measurements = []
        for d in calib_info['data']:
            vis_json, source_json = d
            cv, timestamp = api_imaging.vis_calibrated(vis_json, config, gains, phase_offsets, flag_list)
            src_list = elaz.from_json(source_json, 0.0)
        if not ARGS.show_sources:
            src_list = None
    else:
        logger.info("Getting Data from API: {}".format(ARGS.api))

        api = api_handler.APIhandler(ARGS.api)
        config = api_handler.get_config(api)

        gains = api.get('calibration/gain')

        if (ARGS.vis is None):
            vis_json = api.get('imaging/vis')
        else:
            with open(ARGS.vis, 'r') as json_file:
                vis_json = json.load(json_file)

        ts = api_imaging.vis_json_timestamp(vis_json)
        if ARGS.show_sources:
            source_json = api.get_url(api.catalog_url(config, datestr=ts.isoformat()))

        logger.info("Data Download Complete")

        cv, timestamp = api_imaging.vis_calibrated(vis_json, config, gains['gain'], gains['phase_offset'], flag_list=[])

    api_imaging.rotate_vis(ARGS.rotation, cv, reference_positions = deepcopy(config.get_antenna_positions()))
    
    time_repr = "{:%Y_%m_%d_%H_%M_%S_%Z}".format(timestamp)

    # Processing
    should_make_images = ARGS.display or ARGS.PNG or ARGS.fits
    
    spot = Spotless(cv)
    
    src_list = get_source_list(source_json, el_limit=ARGS.elevation, jy_limit=1e4)
    
    if should_make_images:
        nside = 2**ARGS.nside
        spt = spot.display(plt, src_list, nside, False)
        handle_image(ARGS, None, "" + ARGS.title, time_repr, src_list)
    
    if ARGS.beam:
        spot.beam(plt, nside)
        handle_image(ARGS, None, "dirty", "beam", src_list)
