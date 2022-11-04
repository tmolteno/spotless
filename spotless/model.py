#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#

import numpy as np
import json

from .source import PointSource


class Model(object):

    def __init__(self):
        self.objects = []

    def add_source(self, src):
        self.objects.append(src)

    def __getitem__(self, index):
        return self.objects[index]

    def model_vis(self, u_arr, v_arr, w_arr):
        vis = np.zeros_like(u_arr, dtype=np.complex64)

        for src in self.objects:
            vis += src.get_vis(u_arr, v_arr, w_arr)
        return vis

    def to_dict(self):
        ret = {}
        objs = [obj.to_dict() for obj in self.objects]
        ret["model"] = objs
        return ret

    def __repr__(self):
        return json.dumps(self.to_dict())

    def brightest(self):
        '''Return Brightest Source'''
        ret = None
        b = 0.0
        for src in self.objects:
            if src.a > b:
                b = src.a
                ret = src
        return ret

    def faintest(self):
        '''Return Faintest Source'''
        ret = None
        b = 9e99
        for src in self.objects:
            if src.a < b:
                b = src.a
                ret = src
        return ret

    def to_vector(self):
        ret = []
        for src in self.objects:
            ret.append(src.a)
            ret.append(src.el)
            ret.append(src.az)
        return np.array(ret)

    @classmethod
    def from_vector(cls, vect):
        ret = cls()
        bits = np.split(vect, len(vect)/3)
        for src in bits:
            a, el, az = src
            ret.add_source(PointSource(a, el, az))
        return ret
