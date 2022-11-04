#!/usr/bin/env python
"""
    Calibrate the Telescope from the RESTful API
    using the spotless imaging algorithm

    Tim Molteno 2017.
"""
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
import time

from copy import deepcopy

import itertools

from tart.operation import settings
from tart.imaging import visibility
from tart.imaging import calibration
from tart.imaging import synthesis
from tart.imaging import elaz

from tart.util.angle import from_rad

from tart_tools import api_imaging

from spotless import Spotless, get_source_list

def split_param(x):
    rot_degrees = x[0]
    re = np.concatenate(([1], x[1:24]))
    im = np.concatenate(([0], x[24:47]))
    gains = np.sqrt(re*re + im*im)
    phase_offsets = np.arctan2(im, re)
    return rot_degrees, gains, phase_offsets

def join_param(rot_degrees, gains, phase_offsets):
    ret = np.zeros(47)
    ret[0] = rot_degrees
    z = gains[1:24]*np.exp(phase_offsets[1:24]*1j)
    ret[ 1:24] = z.real
    ret[24:47] = z.imag
    return ret


def param_to_json(x):
    rot_degrees, gains, phase_offsets = split_param(x)
    ret = {'gain':  np.round(gains, 4).tolist(),
           'rot_degrees': rot_degrees,
           'phase_offset': np.round(phase_offsets,4).tolist()}
    return ret

def output_param(x, fp=None):
    ret = param_to_json(x)
    if (fp is None):
        print((json.dumps(ret, indent=4, separators=(',', ': '))))
    else:
        json.dump(ret, fp, indent=4, separators=(',', ': '))


def elaz_to_r(el, az, r=1.0):
    cel = np.cos(el)
    x = cel*np.sin(az)
    y = cel*np.cos(az)
    z = np.sin(el)
    return np.array([x,y,z])*r

def elaz_distance(el_0, az_0, el_1, az_1):
    ''' Return the angle between two positions in
        the horizontal coordinate system
    '''
    assert (el_0 <= np.pi/2)
    assert (el_0 >= -np.pi/2)
    p0 = elaz_to_r(el_0, az_0)
    p1 = elaz_to_r(el_1, az_1)
    d = np.arccos(np.dot(p0, p1))  # |a||b| cos(theta) = a.b
    return d

def calc_score_aux(opt_parameters, measurements, window_deg, original_positions):
    global triplets, ij_index, jk_index, ik_index
    rot_degrees, gains, phase_offsets = split_param(opt_parameters)

    ret_zone = 0.0
    
    ant_idxs = np.arange(24)
    
    for m in measurements:
        cv, ts, src_list = m

        cv.set_phase_offset(ant_idxs, phase_offsets)
        cv.set_gain(ant_idxs, gains)
        api_imaging.rotate_vis(rot_degrees, cv, original_positions)

        spot = Spotless(cv)
        spot.working_nside = 2**5
        spot.deconvolute()
        model_dict = spot.model.to_dict()

        zone_score = 0.0

        for known_src in src_list:
            closest = 6.0
            power = 0.0
            n = 0
            for model_src in spot.model:
                d = elaz_distance(known_src.el_r, known_src.az_r, model_src.el, model_src.az)
                n += 1
                if (d < closest):
                    closest = d
                    #power = model_src.power

            zone_score += 1.0/(closest + np.radians(1))**2

        zone_score = -zone_score/(len(src_list) + n)
        ret_zone += zone_score
    sys.stdout.write("{}\r".format(ret_zone))
    sys.stdout.flush()
    return ret_zone


def load_data_from_json(vis_json, src_json, config, gains, phases, flag_list, el_threshold):

    cv, ts = api_imaging.vis_calibrated(vis_json, config, gains, phases, flag_list)
    src_list = elaz.from_json(src_json, el_threshold)
    return cv, ts, src_list


def calc_score(opt_parameters, config, measurements, window_deg, original_positions, update=False, show=False):
    global N_IT, method, output_directory, f_vs_iteration

    ret = calc_score_aux(opt_parameters, measurements, window_deg, original_positions)

    if N_IT%100==0:
        print(N_IT, ret)
        f_vs_iteration.append(ret)

    N_IT += 1
    return ret

from scipy import optimize
import json

def bfgs_callback(xk):
    print(("Iteration at {}".format(xk)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calibrate the tart telescope from API using spotless.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--api', required=False, default='https://tart.elec.ac.nz/signal', help="Telescope API server URL.")
    parser.add_argument('--file', required=False, default='calibration_data.json', help="Calibration Output JSON file.")
    parser.add_argument('--show', action="store_true", help="show instead of save.")
    parser.add_argument('--get-gains', action="store_true", help="Start from current knowledge of gains.")
    parser.add_argument('--dir', required=False, default='.', help="Output directory.")
    parser.add_argument('--elevation', type=float, default=20.0, help="Elevation threshold for sources]")

    ARGS = parser.parse_args()
    
    # Load calibration data
    with open(ARGS.file, 'r') as json_file:
        calib_info = json.load(json_file)

    info = calib_info['info']
    ant_pos = calib_info['ant_pos']
    config = settings.from_api_json(info['info'], ant_pos)

    flag_list = [] # [4, 5, 14, 22]

    output_directory = ARGS.dir
    f_vs_iteration = []
    
    original_positions = deepcopy(config.get_antenna_positions())

    gains_json = calib_info['gains']
    print(gains_json['gain'])
    gains = np.asarray(gains_json['gain'])
    phase_offsets = np.asarray(gains_json['phase_offset'])
    print(gains)
    print(phase_offsets)
    
    if (ARGS.get_gains):
        api = api_imaging.APIhandler(ARGS.api)
        gains_json = api.get('calibration/gain')
        gains = np.asarray(gains_json['gain'])
        phase_offsets = np.asarray(gains_json['phase_offset'])
    

        
    config = settings.from_api_json(info['info'], ant_pos)
 
    init_parameters = join_param(0.0, gains, phase_offsets)
    output_param(init_parameters)

    measurements = []
    for d in calib_info['data']:
        vis_json, src_json = d
        cv, ts, src_list = load_data_from_json(vis_json, src_json, config, gains, phase_offsets, flag_list, el_threshold=ARGS.elevation)
        measurements.append([cv, ts, src_list])
    
    N_IT = 0
    window_deg = 3.0
    
    s = calc_score(init_parameters, config, measurements, window_deg, original_positions, update=False, show=False)

    f = lambda param: calc_score(param, config, measurements, window_deg, original_positions, update=False, show=False)

    bounds = [(-6, 6)]  # Bounds for the rotation parameter
    for i in range(46):
        bounds.append((-2,2)) # Bounds for all other parameters (real and imaginary components)

    np.random.seed(555)   # Seeded to allow replication.

    bfgs_options = {'ftol':1e-5, 'eps':1e-5}
    ret = optimize.minimize(f, init_parameters, method='L-BFGS-B', bounds=bounds, callback=bfgs_callback, options=bfgs_options)

    rot_degrees = ret.x[0]
    output_json = param_to_json(ret.x)
    output_json['message'] = ret.message
    output_json['optimum'] = ret.fun
    output_json['iterations'] = ret.nit

    new_positions =  settings.rotate_location(rot_degrees, np.array(original_positions).T)
    print(new_positions)
    pos_list = (np.array(new_positions).T).tolist()
    print(pos_list)
    output_json['antenna_positions'] = pos_list
    
    with open('{}/{}_opt_json.json'.format(output_directory, 'spotless'), 'w') as fp:
        json.dump(output_json, fp, indent=4, separators=(',', ': '))

    f_history_json = {}
    f_history_json['start'] = s
    f_history_json['history'] = f_vs_iteration

    with open('{}/{}_history.json'.format(output_directory, 'spotless'), 'w') as fp:
        json.dump(f_history_json, fp, indent=4, separators=(',', ': '))
