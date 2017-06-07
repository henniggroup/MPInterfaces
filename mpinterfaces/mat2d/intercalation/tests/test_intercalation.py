import os
import unittest

from pymatgen import Structure

from mpinterfaces.mat2d.intercalation.startup import *
from mpinterfaces.mat2d.intercalation.analysis import *


__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "tests"))


class AnalysisTest(unittest.TestCase):

    def test_get_interstitial_sites(self):
        os.chdir(ROOT)
        os.chdir("MoS2")
        structure = Structure.from_file("POSCAR")
        test_interstitial = get_interstitial_sites(structure)["tetrahedral"][0]
        control_interstitial = ((3.29906811,-0.77121991,17.80937681), 4,
                                1.34201354)
        for i in range(3):
            self.assertAlmostEqual(test_interstitial[0][i],
                                   control_interstitial[0][i])


class StartupTest(unittest.TestCase):

    def test_inject_ions(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        structure = Structure.from_file('POSCAR')
        structure = inject_ions(structure, 'Li', 0.25)
        structure.to(fmt='POSCAR', filename='test_intercalated_POSCAR')
        control_lines = open('control_intercalated_POSCAR').readlines()
        test_lines = open('test_intercalated_POSCAR').readlines()
        for i in range(len(control_lines)):
            self.assertEqual(control_lines[i], test_lines[i])

if __name__ == '__main__':
    unittest.main()
