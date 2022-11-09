#
# Copyright Tim Molteno 2017-2022 tim@elec.ac.nz
# License GPLv3
#

import unittest

import numpy as np

from spotless import sphere


class TestUtil(unittest.TestCase):

    def setUp(self):
        pass

    # def test_hp_elaz_hp(self):
    #     theta = np.random.rand(100)*np.pi/2
    #     phi = np.random.rand(100)*np.pi*2
    # 
    #     el, az = sphere.HealpixSphere.hp2elaz(theta, phi)
    # 
    #     theta2, phi2 = sphere.HealpixSphere.elaz2hp(el, az)
    # 
    #     for i in range(100):
    #         self.assertAlmostEqual(theta2[i], theta[i])
    #         self.assertAlmostEqual(phi2[i], phi[i])
    # 
    # def test_elaz(self):
    #     el = np.pi/2
    #     az = 0.0
    # 
    #     theta, phi = sphere.HealpixSphere.elaz2hp(el, az)
    #     self.assertEqual(theta, 0)
    #     self.assertEqual(phi, 0)
