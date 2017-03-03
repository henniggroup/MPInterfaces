import unittest

import os

from pymatgen.core.structure import Structure

from twod_materials.intercalation.startup import *
from twod_materials.intercalation.analysis import *

import twod_materials

__author__ = "Michael V. Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__version__ = "1.6"
__maintainer__ = "Michael V. Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')
ROOT = os.path.join(PACKAGE_PATH, 'intercalation/tests')

class StartupTest(unittest.TestCase):

    """ This will be hard to test, since inject_ions requires Zeo++
        to be installed...

    def test_inject_ions(self):
        os.chdir(ROOT)
        os.chdir('MoS2')
        structure = Structure.from_file('POSCAR')
        structure = inject_ions(structure, 'Li', 0.25)
        structure.to(fmt='POSCAR', filename='test_intercalated_POSCAR')
        control_lines = open('control_intercalated_POSCAR').readlines()
        test_lines = open('control_intercalated_POSCAR').readlines()
        for i in range(len(control_lines)):
            self.assertEqual(control_lines[i], test_lines[i])
    """

if __name__ == '__main__':
    unittest.main()
