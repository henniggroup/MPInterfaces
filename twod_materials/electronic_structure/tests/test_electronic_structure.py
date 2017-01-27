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

    def test_get_band_edges(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        test_edges = get_band_edges()
        control_edges = {u'dn_cbm': -4.5298536713404509,
                         u'dn_vbm': -6.2530536713404512,
                         u'efermi': -6.2162497913404513,
                         u'up_cbm': -4.5298536713404509,
                         u'up_vbm': -6.2530536713404512}
        self.assertEqual(test_edges, control_edges)

    """ These tests are annoying because they require LaTeX, which
        has to be installed on Travis and takes forever.

    def test_plot_band_alignments_creates_file(self):
        os.chdir(ROOT)
        plot_band_alignments(['MoS2'])
        self.assertTrue(os.path.isfile('band_alignments.pdf'))
        os.system('rm band_alignments.pdf')


    def test_plot_local_potential_creates_file(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        plot_local_potential()
        self.assertTrue(os.path.isfile('locpot.pdf'))
        os.system('rm locpot.pdf')


    def test_plot_band_structure_creates_file(self):
        os.chdir(ROOT)
        os.chdir('band_structure_control')
        plot_band_structure()
        self.assertTrue(os.path.isfile('band_structure.pdf'))
        os.system('rm band_structure.pdf')


    def test_plot_color_projected_bands_creates_file(self):
        os.chdir(ROOT)
        os.chdir('band_structure_control')
        plot_color_projected_bands()
        self.assertTrue(os.path.isfile('color_projected_bands.pdf'))
        os.system('rm color_projected_bands.pdf')


    def test_plot_elt_projected_bands_creates_file(self):
        os.chdir(ROOT)
        os.chdir('band_structure_control')
        plot_elt_projected_bands()
        self.assertTrue(os.path.isfile('elt_projected_bands.pdf'))
        os.system('rm elt_projected_bands.pdf')


    def test_plot_orb_projected_bands_creates_file(self):
        os.chdir(ROOT)
        os.chdir('band_structure_control')
        orbitals = {'B': ['s', 'p'], 'N': ['s', 'p']}
        plot_orb_projected_bands(orbitals)
        self.assertTrue(os.path.isfile('orb_projected_bands.pdf'))
        os.system('rm orb_projected_bands.pdf')
    """

if __name__ == '__main__':
    unittest.main()
