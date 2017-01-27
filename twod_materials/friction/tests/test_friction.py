import unittest

import os

from twod_materials.friction.startup import *
from twod_materials.friction.analysis import *

import twod_materials


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')
ROOT = os.path.join(PACKAGE_PATH, 'friction/tests')

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
        control_data = {'F_N': [10.917128940360275],
                        'mu': [0.71346395371219118],
                        'F_f': [7.7889779769752261]}
        for key in control_data:
            self.assertEqual(test_data[key][0], control_data[key][0])


if __name__ == '__main__':
    unittest.main()
