#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#

import disko
import logging

import numpy as np
import json
from . import sphere

from tart.util import constants

logger = logging.getLogger(__name__)

class PointSource(object):
    '''
    A single point source in el,az co-ordinates
    '''

    def __init__(self, a, el, az):
        self.a = a
        self.el = el
        self.az = az
        self.power = None   # The power will be set after fitting

    def get_power(self):
        if self.power is None:
            return self.a
        return self.power

    def to_dict(self):
        ret = {}
        ret["a"] = self.a
        ret["el"] = np.degrees(self.el)
        ret["az"] = np.degrees(self.az)
        if (self.power is not None):
            ret["p"] = float(self.power)
            # return "{p:None, a:{:03.2f}, el:{:03.1f}, az:{:03.1f}},".format(self.a, np.degrees(self.el), np.degrees(self.az))
        # else:
            # return "{p:{:03.2f}, a:{:03.2f}, el:{:03.1f}, az:{:03.1f}},".format(self.power,
            # self.a, np.degrees(self.el), np.degrees(self.az))
        return ret

    def __repr__(self):
        return json.dumps(self.to_dict())

    def get_vis(self, u_arr, v_arr, w_arr):
        '''
            Generate a list of visibilities from the source at u,v,w points
            See "The non-coplanar baselines effect in radio interferometry: The W-projection algorithm"
            by Tim Cornwell. Note that this paper chooses direction cosines (l,m) as its basis for the
            integral and the image coordinates. It introduces a factor sqrt(1 - l**2 - m**2) which
            implies that visibilities are infinite when the object is on the horizon (l or m ==1). This
            is really just a correction because dl dm approaches zero at this point.

            The Smirnov papers on the RIME provide a formalism that is easier to derive algorithms for.

            TODO: Use an antenna model (for the whole telescope initially, and then eventually for each antenna)
            to modify the source amplitude (self.a) and predict the visibilities from a low-elevation source.
        '''
        l, m, n = sphere.elaz2lmn(self.el, self.az)

        p2j = disko.jomega(constants.L1_FREQ)
        vis = self.a*self.a * \
            np.exp(-p2j*(u_arr*l + v_arr*m + w_arr*(n - 1.0)))
        return vis

    def get_bounds(self, d_el, el_threshold_r=0):
        '''
            TODO check that the source remains inside the sphere.
        '''
        logger.info(f" get_bounds({d_el}, {el_threshold_r})")
        d_az = np.abs(d_el/(np.cos(self.el) + 0.001))
        
        el_lower = max(self.el - d_el, el_threshold_r)
        el_upper = max(el_threshold_r, self.el + d_el)
        
        logger.info(f" Elevation: ({el_lower}, {el_upper})")
        
        return [(0.0, None),
                (el_lower, min(el_upper, np.pi/2)),
                (self.az - d_az, self.az + d_az)]
