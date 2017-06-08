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

    # In practice it is very difficult to reproduce the
    # number of interstitial sites and their exact
    # locations on different architectures. This is not
    # a big problem and is the reason the tests are a little
    # vague.
    def test_get_interstitial_sites_with_octahedra(self):
        os.chdir(ROOT)
        structure = Structure.from_file("POSCAR_Cu")
        test_ints = get_interstitial_sites(structure, octahedra=True)
        self.assertTrue(len(test_ints["tetrahedral"]) != 0)
        self.assertTrue(len(test_ints["hexahedral"]) != 0)
        self.assertTrue(len(test_ints["octahedral"]) != 0)


    def test_get_interstitial_sites_without_octahedra(self):
        os.chdir(ROOT)
        structure = Structure.from_file("POSCAR_Cu")
        test_ints = get_interstitial_sites(structure, octahedra=False)
        self.assertTrue(len(test_ints["tetrahedral"]) != 0)


class StartupTest(unittest.TestCase):

    def test_inject_ions(self):
        os.chdir(ROOT)
        structure = Structure.from_file('POSCAR_MoS2')
        structure = inject_ions(structure, 'Li', 0.25)
        self.assertTrue(structure.num_sites == 4)

if __name__ == '__main__':
    unittest.main()
