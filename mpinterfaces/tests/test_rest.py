from __future__ import unicode_literals

import unittest

import os

import json

from mpinterfaces.rest import MWRester
from pymatgen.core.structure import Structure

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


TEST_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "test_files")
MWR = MWRester()

"""
class RestTest(unittest.TestCase):

    def test_get_data_for_mw_43(self):
        test_data = MWR.get_data('mw-43')[0]
        with open(os.path.join(TEST_FILES_DIR, 'mw-43_data.json'), 'r') as j:
            control_data = json.load(j)['response'][0]
        for k in ['last_updated', 'unit_cell_formula', 'material_id']:
            self.assertEqual(test_data[k], control_data[k])

    def test_get_structure_by_material_id_for_mw_43(self):
        test_structure = MWR.get_structure_by_material_id('mw-43')
        control_structure = Structure.from_file(
            os.path.join(TEST_FILES_DIR, 'POSCAR_mw-43')
        )
        self.assertAlmostEqual(test_structure.lattice.a,
                               control_structure.lattice.a, places=3)
        self.assertAlmostEqual(test_structure.lattice.b,
                               control_structure.lattice.b, places=3)
        self.assertAlmostEqual(test_structure.lattice.c,
                               control_structure.lattice.c, places=3)
        self.assertAlmostEqual(test_structure.sites[0],
                               control_structure.sites[0], places=3)


if __name__ == '__main__':
    unittest.main()
"""
