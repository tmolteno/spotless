#
# Copyright Tim Molteno 2017-2026 tim@elec.ac.nz
# License GPLv3
#

import logging

import numpy as np
from scipy.optimize import minimize

from .model import Model
from .source import PointSource
from .spotless import Spotless

logger = logging.getLogger(__name__)


class MultiSpotless(Spotless):
    def __init__(self, disko, sphere):
        super(MultiSpotless, self).__init__(disko, sphere)

    def step(self):
        """
        Return a list of residual visibilities and a source location
        so that the sum of the residual and the source visibilities
        is conserved.
        """
        d_el = self.disko.get_beam_width().radians() / 2

        m_vis = self.model.model_vis(
            self.disko.u_arr, self.disko.v_arr, self.disko.w_arr
        )

        a_0, el_0, az_0, p0 = self.estimate_initial_point_source(self.vis_arr - m_vis)

        a_init = max(a_0, 0.01)
        model_vect = self.model.to_vector()
        x0 = np.append(model_vect, [a_init, el_0, az_0])

        # Get the bounds of the existing points. Each dimension is
        # constrained appropriately.
        bounds = []
        for src in self.model:
            for b in src.get_bounds(d_el, self.el_threshold_r):
                bounds.append(b)

        src = PointSource(a_init, el_0, az_0)
        for b in src.get_bounds(d_el, self.el_threshold_r):
            bounds.append(b)

        # Use Nelder-Mead for MultiSpotless: the mixed amplitude/angle
        # parameter space is poorly scaled, making gradient-based
        # methods (L-BFGS-B) unreliable. Nelder-Mead only uses function
        # values and handles this naturally.
        #
        # Set generous maxfev/maxiter since dimensionality grows with
        # source count.  If the optimizer fails, retry once with a
        # perturbed starting point.
        n_params = len(x0.flatten())
        maxfev = max(n_params * 500, 5000)
        fmin = minimize(
            self.f_n,
            x0.flatten(),
            method="Nelder-Mead",
            bounds=bounds,
            options={
                "xatol": 1e-3,
                "fatol": 0.1,
                "adaptive": True,
                "maxfev": maxfev,
                "maxiter": maxfev,
            },
        )

        # Retry with a perturbed starting point if the optimizer
        # failed to converge (simplex hit maxfev or collapsed).
        if not fmin.success:
            logger.warning(
                "Nelder-Mead did not converge (nfev=%d), retrying...",
                fmin.nfev,
            )
            x0_retry = x0.flatten() * (1.0 + 0.1 * np.random.randn(len(x0.flatten())))
            x0_retry = np.clip(
                x0_retry,
                [b[0] if b[0] is not None else -np.inf for b in bounds],
                [b[1] if b[1] is not None else np.inf for b in bounds],
            )
            fmin = minimize(
                self.f_n,
                x0_retry,
                method="Nelder-Mead",
                bounds=bounds,
                options={
                    "xatol": 1e-3,
                    "fatol": 0.1,
                    "adaptive": True,
                    "maxfev": maxfev,
                    "maxiter": maxfev,
                },
            )
        self._opt_nfev = fmin.nfev
        self._opt_nit = fmin.nit
        self._opt_success = fmin.success
        residual_power = fmin.fun

        self.model = Model.from_vector(fmin.x)

        self.residual_vis = self.vis_arr - self.model.model_vis(
            self.disko.u_arr, self.disko.v_arr, self.disko.w_arr
        )

        return self.model, residual_power, p0

    def f_n(self, x):
        """
        Find the power in the residual
        """
        mod = Model.from_vector(x)
        m_vis = mod.model_vis(self.disko.u_arr, self.disko.v_arr, self.disko.w_arr)
        return self.power(self.vis_arr - m_vis)
