#
# Copyright Tim Molteno 2017 tim@elec.ac.nz
#

import unittest

import numpy as np

from spotless import source
from spotless import sphere


class TestSource(unittest.TestCase):

    def setUp(self):
        pass

    def test_vis(self):
        src = source.PointSource(1.0, np.pi/2, 0.0)

        # TODO
