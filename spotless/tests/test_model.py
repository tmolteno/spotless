#
# Copyright Tim Molteno 2017 tim@elec.ac.nz
#


import unittest

from spotless import Model
from spotless import PointSource


class TestModel(unittest.TestCase):

    def setUp(self):
        self.spot = Model()
        self.spot.add_source(PointSource(0.1, 0.11, 0.12))
        self.spot.add_source(PointSource(0.2, 0.3, 0.4))

    def test_brightest(self):
        bright = self.spot.brightest()
        self.assertEqual(bright.a, 0.2)

    def test_faintest(self):
        bright = self.spot.faintest()
        self.assertEqual(bright.a, 0.1)

    def test_to_vector(self):
        vect = self.spot.to_vector()
        self.assertEqual(vect[0], 0.1)
        self.assertEqual(vect[-1], 0.4)
        self.assertEqual(vect.shape, (6,))
        
    def test_to_from_vector(self):
        vect = self.spot.to_vector()
        spot2 = Model.from_vector(vect)
        for a, b in zip(vect, spot2.to_vector()):
            self.assertEqual(a, b)
