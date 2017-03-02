import unittest

import os

from twod_materials.electronic_structure.startup import *
from twod_materials.electronic_structure.analysis import *

import twod_materials


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')
ROOT = os.path.join(PACKAGE_PATH, 'electronic_structure/tests')


class StartupTest(unittest.TestCase):

    def test_run_pbe_calculation_creates_files(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        run_pbe_calculation(submit=False)
        for f in ['POSCAR', 'INCAR', 'KPOINTS']:
            self.assertTrue(os.path.isfile('pbe_bands/{}'.format(f)))
        os.system('rm -r pbe_bands')


    def test_run_hse_calculation_creates_files(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        run_hse_calculation(submit=False)
        for f in ['POSCAR', 'INCAR', 'KPOINTS', 'runjob']:
            self.assertTrue(os.path.isfile('hse_bands/{}'.format(f)))
        os.system('rm -r hse_bands')


    def test_run_hse_prep_calculation_creates_files(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        run_hse_prep_calculation(submit=False)
        for f in ['POSCAR', 'INCAR', 'KPOINTS', 'runjob']:
            self.assertTrue(os.path.isfile('hse_prep/{}'.format(f)))
        os.system('rm -r hse_prep')


class AnalysisTest(unittest.TestCase):

    def test_get_band_edges_for_MoS2(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        test_edges = get_band_edges()
        control_edges = {u'dn_cbm': -4.5298536713404509,
                         u'dn_vbm': -6.2530536713404512,
                         u'efermi': -6.2162497913404513,
                         u'up_cbm': -4.5298536713404509,
                         u'up_vbm': -6.2530536713404512}
        for edge in test_edges:
            self.assertAlmostEqual(test_edges[edge], control_edges[edge])


    def test_get_fermi_velocities_for_MoS2_and_FeCl2(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        test_velocities = get_fermi_velocities()
        control_velocities = []
        self.assertEqual(test_velocities, control_velocities)
        # Now FeCl2, which does have metallic bands.
        os.chdir('../FeCl2')
        test_velocities = get_fermi_velocities()
        control_velocities = [1211920261158078.5, 1213975748043003.5,
                              3440582787448865.5, 1525493730297101.2,
                              2953560839174396.5, 2212830757399824.8,
                              2855773386125549.0, 1411440366276435.5,
                              2970869118108168.0, 2018812598303622.8]
        for i in range(len(test_velocities)):
            self.assertEqual(test_velocities[i], control_velocities[i])


    def test_plot_band_alignments_creates_data(self):
        os.chdir(ROOT)
        ax = plot_band_alignments(['MoS2'], fmt='None')
        self.assertEqual(ax.get_children()[0].get_height(), 1.23)


    def test_plot_local_potential_creates_data(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        ax = plot_local_potential(fmt='None')
        test_line = ax.get_lines()[0]
        test_x, test_y = test_line.get_xdata(), test_line.get_ydata()
        control_lines = open('locpot_data.txt').readlines()
        control_x = [float(d) for d in control_lines[0].split(',')]
        control_y = [float(d) for d in control_lines[1].split(',')]
        for i in range(len(control_x)):
            self.assertAlmostEqual(test_x[i], control_x[i])
            self.assertAlmostEqual(test_y[i], control_y[i])


    def test_plot_band_structure_creates_data(self):
        os.chdir(ROOT)
        os.chdir('band_structure_control')
        ax = plot_band_structure(fmt='None')
        test_line = ax.get_lines()[0]
        test_x, test_y = test_line.get_xdata(), test_line.get_ydata()
        control_x = [
            0.0, 0.07230748, 0.14461496, 0.21692245, 0.28922993, 0.36153741,
            0.43384489, 0.50615237, 0.57845986, 0.65076734, 0.72307482,
            0.7953823 , 0.86768978, 0.93999727, 1.01230475, 1.08461223,
            1.15691971, 1.22922719, 1.30153468, 1.37384216, 1.44614964
        ]
        control_y = [
            -19.4417, -19.4287, -19.3895, -19.3245, -19.2338, -19.1181,
            -18.9778, -18.8139, -18.6275, -18.4202, -18.1937, -17.9507,
            -17.6945, -17.4296, -17.162 , -16.8998, -16.6534, -16.4361,
            -16.2634, -16.1514, -16.1125
        ]
        for i in range(len(control_x)):
            self.assertAlmostEqual(test_x[i], control_x[i])
            self.assertAlmostEqual(test_y[i], control_y[i])


    def test_plot_color_projected_bands_creates_file(self):
        os.chdir(ROOT)
        os.chdir('band_structure_control')
        ax = plot_color_projected_bands(fmt='None')
        test_line = ax.get_lines()[6]
        test_x, test_y = test_line.get_xdata(), test_line.get_ydata()
        control_x, control_y = [0.14461496, 0.21692245], [-19.3895, -19.3245]
        for i in range(len(control_x)):
            self.assertAlmostEqual(test_x[i], control_x[i])
            self.assertAlmostEqual(test_y[i], control_y[i])


    def test_plot_elt_projected_bands_creates_data(self):
        os.chdir(ROOT)
        os.chdir('band_structure_control')
        ax = plot_elt_projected_bands(fmt='None')
        test_line = ax.get_lines()[7]
        test_x, test_y = test_line.get_xdata()[0], test_line.get_ydata()[0]
        control_x, control_y = 0.14461496, -19.3895
        self.assertAlmostEqual(test_x, control_x)
        self.assertAlmostEqual(test_y, control_y)


    def test_plot_orb_projected_bands_creates_data(self):
        os.chdir(ROOT)
        os.chdir('band_structure_control')
        orbitals = {'B': ['s', 'p'], 'N': ['s', 'p']}
        ax = plot_orb_projected_bands(orbitals, fmt='None')
        test_line = ax.get_lines()[8]
        test_x, test_y = test_line.get_xdata()[0], test_line.get_ydata()[0]
        control_x, control_y = 0.21692245, -19.3245
        self.assertAlmostEqual(test_x, control_x)
        self.assertAlmostEqual(test_y, control_y)


if __name__ == '__main__':
    unittest.main()
