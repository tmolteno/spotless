#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#

import numpy as np
import logging

from scipy.optimize import minimize

from .spotless import Spotless
from .model import Model
from .source import PointSource

logger = logging.getLogger(__name__)


class MultiSpotless(Spotless):

    def __init__(self, disko, sphere):
        super(MultiSpotless, self).__init__(disko, sphere)

    def step(self):
        '''
            Return a list of residual visibilities and a source location
            so that the sum of the residual and the source visibilities
            is conserved.
        '''
        d_el = self.disko.get_beam_width().radians()/2

        m_vis = self.model.model_vis(self.disko.u_arr,
                                     self.disko.v_arr, self.disko.w_arr)

        a_0, el_0, az_0, p0 = self.estimate_initial_point_source(
            self.vis_arr - m_vis)

        model_vect = self.model.to_vector()
        x0 = np.append(model_vect, [0.1, el_0, az_0])

        # Get the bounds of the existing points. Each dimension is constrained appropriately
        bounds = []
        for src in self.model:
            for b in src.get_bounds(d_el, self.el_threshold_r):
                bounds.append(b)

        src = PointSource(0.1, el_0, az_0)
        for b in src.get_bounds(d_el, self.el_threshold_r):
            bounds.append(b)

        # fmin = minimize(self.f_n, x0.flatten(),
        #                 method='L-BFGS-B', bounds=bounds)
        fmin = minimize(self.f_n, x0.flatten(),
                        method='Nelder-mead', bounds=bounds)
        residual_power = fmin.fun

        self.model = Model.from_vector(fmin.x)

        self.residual_vis = self.vis_arr - \
            self.model.model_vis(self.disko.u_arr, self.disko.v_arr, self.disko.w_arr)

        return self.model, residual_power, p0

    def f_n(self, x):
        '''
            Find the power in the residual
        '''
        mod = Model.from_vector(x)
        m_vis = mod.model_vis(self.disko.u_arr,
                              self.disko.v_arr,
                              self.disko.w_arr)
        return self.power(self.vis_arr - m_vis)
