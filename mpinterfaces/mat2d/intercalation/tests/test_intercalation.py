import os
import unittest

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "tests"))


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
    pass

if __name__ == '__main__':
    unittest.main()
