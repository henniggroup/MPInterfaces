import os
import unittest

from mpinterfaces import MPINT_CONFIG
from mpinterfaces.mat2d.friction import *
from mpinterfaces.mat2d.friction.analysis import get_basin_and_peak_locations, \
    get_mu_vs_F_N

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "tests"))


class StartupTest(unittest.TestCase):

    def test_run_gamma_calculations(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        run_gamma_calculations(submit=False)
        self.assertTrue(os.path.isfile('friction/lateral/0x0/POSCAR'))
        os.system('rm -r friction')

    def test_run_normal_force_calculations(self):
        os.chdir(ROOT)
        os.chdir('MoS2_with_lateral')
        run_normal_force_calculations(('0x0', '2x5'), submit=False)
        self.assertTrue(os.path.isfile('friction/normal/1.5/0x0/POSCAR'))
        os.system('rm -r friction/normal')


class AnalysisTest(unittest.TestCase):

    def test_get_basin_and_peak_locations(self):
        os.chdir(ROOT)
        os.chdir('MoS2_with_lateral')
        self.assertEqual(get_basin_and_peak_locations(), (('2x5', '0x0')))

    def test_get_mu_vs_F_N(self):
        os.chdir(ROOT)
        os.chdir('MoS2_with_lateral_and_normal')
        test_data = get_mu_vs_F_N('2x5')
        control_data = {'F_N': [2.7292822350900687],
                        'mu': [0.71346395371219118],
                        'F_f': [1.9472444942438065]}
        for key in control_data:
            self.assertEqual(test_data[key][0], control_data[key][0])


if __name__ == '__main__':
    unittest.main()
