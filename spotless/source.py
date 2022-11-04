#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#

import disko

import numpy as np
import json
from . import sphere

from tart.util import constants

class PointSource(object):
    '''
    A single point source in el,az co-ordinates
    '''

    def __init__(self, a, el, az):
        self.a = a
        self.el = el
        self.az = az
        self.power = None   # The power will be set after fitting

    def to_dict(self):
        ret = {}
        ret["a"] = self.a
        ret["el"] = np.degrees(self.el)
        ret["az"] = np.degrees(self.az)
        if (self.power is not None):
            ret["p"] = self.power
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

    def get_bounds(self, d_el):
        d_az = d_el/(np.cos(self.el) + 0.001)
        return [(0.0, 1.0),
                (max(self.el - d_el, 0), min(self.el + d_el, np.pi/2)),
                (self.az - d_az, self.az + d_az)]
