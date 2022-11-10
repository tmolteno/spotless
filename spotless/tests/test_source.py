#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#

import unittest

import numpy as np

from spotless import source


class TestSource(unittest.TestCase):

    def setUp(self):
        pass

    def test_vis(self):
        src = source.PointSource(1.0, np.pi/2, 0.0)

        # TODO
