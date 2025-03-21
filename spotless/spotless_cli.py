#!/usr/bin/env python
#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#

import matplotlib
import os
if os.name == 'posix' and "DISPLAY" not in os.environ:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
import json
import logging
import sys

from importlib.metadata import version
from copy import deepcopy
from disko import DiSkO, sphere_from_args, sphere_args_parser, disko_from_ms

import numpy as np

from tart.operation import settings

from tart_tools import api_handler
from tart_tools import api_imaging
from tart.imaging import elaz

from .spotless import Spotless
from .spotless import get_source_list
from .multi_spotless import MultiSpotless

from tart2ms import get_array_location

logger = logging.getLogger()


def handle_image(args, img, title, time_repr, src_list=None, sphere=None):
    """ This function manages the output of an image, drawing sources e.t.c."""
    image_title = f"{args.title}_{title}_{time_repr}"
    plt.title(image_title)
    if args.fits:
        fname = '{}.fits'.format(image_title)
        fpath = os.path.join(args.dir, fname)
        api_imaging.save_fits_image(img, fname=fname, out_dir=args.dir, timestamp=time_repr)
        print("Generating {}".format(fname))
    if args.PNG:
        fname = '{}.png'.format(image_title)
        fpath = os.path.join(args.dir, fname)
        plt.savefig(fpath, dpi=300)
        print("Generating {}".format(fname))
    if args.SVG:
        fname = f"{image_title}.svg"
        fpath = os.path.join(args.dir, fname)
        sphere.to_svg(fname=fname, show_grid=True, src_list=src_list, title=image_title)
        print("Generating {}".format(fname))
    if args.PDF:
        fname = '{}.pdf'.format(image_title)
        fpath = os.path.join(args.dir, fname)
        plt.savefig(fpath, dpi=600)
        print("Generating {}".format(fname))
    if args.display:
        plt.show()


def main():
    sphere_parsers = sphere_args_parser()

    parser = argparse.ArgumentParser(
        description='Generate a SPOTLESS image for a TART radio telescope.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=sphere_parsers)

    parser.add_argument('--api', required=False,
                        default='https://tart.elec.ac.nz/signal',
                        help="Telescope API server URL.")
    parser.add_argument('--catalog', required=False,
                        default='https://tart.elec.ac.nz/catalog',
                        help="Catalog API URL.")

    data_group = parser.add_mutually_exclusive_group()
    data_group.add_argument('--file', required=False,
                            default=None,
                            help="Snapshot observation saved JSON file.")

    data_group.add_argument('--ms', required=False,
                            default=None, help="visibility file")
    data_group.add_argument(
        '--vis', required=False, default=None,
        help="Use a local JSON file containing the visibilities.")

    parser.add_argument('--nvis', type=int, default=1000,
                        help="Number of visibilities to use.")
    parser.add_argument('--channel', type=int, default=0,
                        help="Use this frequency channel.")
    parser.add_argument('--field', type=int, default=0,
                        help="Use this FIELD_ID from the measurement set.")
    parser.add_argument('--ddid', type=int, default=0,
                        help="Use this DDID from the measurement set.")

    parser.add_argument('--dir', required=False, default='.',
                        help="Output directory.")

    parser.add_argument('--multimodel', action="store_true",
                        help="Use the SPOTLESS algorithm with multi-dimensional model.")

    parser.add_argument('--beam', action="store_true",
                        help="Generate a gridless beam.")

    parser.add_argument('--display', action="store_true",
                        help="Display Image to the user")
    parser.add_argument('--fits', action="store_true",
                        help="Generate a FITS format image")
    parser.add_argument('--PNG', action="store_true",
                        help="Generate a PNG format image")
    parser.add_argument('--HDF', required=False, default=None,
                        help="Generate an HDF format representation of the field of view")
    parser.add_argument('--PDF', action="store_true",
                        help="Generate a PDF format image")
    parser.add_argument('--SVG', action="store_true",
                        help="Generate a SVG format image")
    parser.add_argument('--show-sources',
                        action="store_true",
                        help="Show known sources on images (only works on PNG).")
    parser.add_argument('--elevation', type=float, 
                        default=20.0, 
                        help="Elevation limit for displaying sources (degrees)")
    parser.add_argument('--title', required=False, 
                        default="spotless", help="Prefix the output files.")

    parser.add_argument('--version', action="store_true",
                        help="Display the current version")
    parser.add_argument('--debug', action="store_true",
                        help="Display debugging information")

    parser.add_argument('--show-model', action="store_true", help="Show the location of the model sources.")

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
        print(f"spotless: Version {version}")
        print("          (c) 2022-2023 Tim Molteno")
        sys.exit(0)

    sphere = sphere_from_args(ARGS)

    if ARGS.file:
        print("Getting Data from file: {}".format(ARGS.file))
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
        disko = DiSkO.from_cal_vis(cv)
        lat = config.get_lat()
        lon = config.get_lon()
        height = config.get_alt()

    elif ARGS.ms:
        print(f"Getting Data from MS file: {ARGS.ms} to {sphere}")

        if not os.path.exists(ARGS.ms):
            raise RuntimeError(f"Measurement set {ARGS.ms} not found")

        min_res = sphere.min_res()
        logger.info(f"Min Res {min_res}")
        disko = disko_from_ms(ARGS.ms, ARGS.nvis, res=min_res,
                              channel=ARGS.channel,
                              field_id=ARGS.field,
                              ddid=ARGS.ddid)

        timestamp = disko.timestamp
# 
        json_info = get_array_location(ARGS.ms)
        lat = json_info['lat']
        lon = json_info['lon']
        height = json_info['height']

    else:
        print("Getting Data from API: {}".format(ARGS.api))

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
            cat_url = api.catalog_url(lon=config.get_lon(),
                                      lat=config.get_lat(),
                                      datestr=ts.isoformat())
            source_json = api.get_url(cat_url)

        print("Data Download Complete")

        cv, timestamp = api_imaging.vis_calibrated(vis_json, config,
                                                   gains['gain'], 
                                                   gains['phase_offset'],
                                                   flag_list=[])
        disko = DiSkO.from_cal_vis(cv)

        lat = config.get_lat()
        lon = config.get_lon()
        height = config.get_alt()

    # CASAcore UVW is conjugated, so to make things consistent with data
    # streaming off telescope we need the vis flipped about
    if ARGS.ms:
        pass # disko.vis_arr = disko.vis_arr.conjugate()
    elif ARGS.file:
        disko.vis_arr = disko.vis_arr.conjugate()
    else:
        pass

    sphere.set_info(timestamp=timestamp,
                    lon=lon, lat=lat, height=height)

    time_repr = "{:%Y_%m_%d_%H_%M_%S_%Z}".format(timestamp)

    # Processing
    should_make_images = ARGS.display or ARGS.PNG or ARGS.fits or ARGS.PDF or ARGS.SVG or ARGS.HDF

    if ARGS.multimodel:
        spot = MultiSpotless(disko, sphere)
    else:
        spot = Spotless(disko, sphere)

    src_list = get_source_list(source_json, el_limit=ARGS.elevation, jy_limit=1e4)

    if should_make_images:
        plt.figure(figsize=(6, 6))
        spt = spot.display(plt=plt, src_list=src_list, sphere=sphere, show_model=ARGS.show_model)
        handle_image(ARGS, None, "gridless",
                     time_repr, src_list, sphere)

    if ARGS.beam:
        spot.beam(plt, nside)
        handle_image(ARGS, None, "beam",
                     time_repr, src_list, sphere)

    spot.deconvolute()

    if should_make_images:
        residual_sphere = sphere.copy()
        spot.image_visibilities(spot.residual_vis, residual_sphere)
        spot.display(plt, src_list, sphere, ARGS.show_model)
        handle_image(ARGS, None, "residual", time_repr, src_list, residual_sphere)

    reconstructed_sphere, src_power, residual_power = spot.reconstruct()
    if should_make_images:
        spot.plot(plt, reconstructed_sphere, 
                  src_list, ARGS.show_model)
        handle_image(ARGS, None, "spotless", time_repr,
                     src_list, reconstructed_sphere)
    if ARGS.HDF:
        hdf_path = os.path.join(ARGS.dir, ARGS.HDF)
        reconstructed_sphere.to_hdf(hdf_path)

    reconstructed_sphere, src_power, residual_power = spot.reconstruct_direct()
    if should_make_images:
        spot.plot(plt, reconstructed_sphere, src_list,
                  ARGS.show_model)
        handle_image(ARGS, None, "spotless_direct",
                     time_repr, src_list,
                     reconstructed_sphere)

    model_dict = spot.model.to_dict()
    model_dict["time"] = timestamp.isoformat()
    if ARGS.display:
        text = json.dumps(model_dict, sort_keys=True, indent=4,
                          ensure_ascii=True)
        logger.info("Model {}".format(text))

    model_name = f"{ARGS.title}_model_{time_repr}.json"
    fpath = os.path.join(ARGS.dir, model_name)

    with open(fpath, 'w') as outfile:
        json.dump(model_dict, outfile,
                  sort_keys=True, indent=4,
                  ensure_ascii=True)
