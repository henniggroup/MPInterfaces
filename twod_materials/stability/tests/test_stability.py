import unittest

import os

from twod_materials.stability.startup import *
from twod_materials.stability.analysis import *

import twod_materials


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')
ROOT = os.path.join(PACKAGE_PATH, 'stability/tests')

class StartupTest(unittest.TestCase):

    def test_relax_creates_files(self):
        os.chdir(ROOT)
        os.chdir('BiTeCl')
        relax(submit=False, force_overwrite=True)
        for f in ['KPOINTS', 'INCAR', 'runjob']:
            self.assertTrue(os.path.isfile(f))
            os.system('rm {}'.format(f))


class AnalysisTest(unittest.TestCase):

    def test_get_hull_distance_for_BiTeCl(self):
        os.chdir(ROOT)
        os.chdir('BiTeCl')
        self.assertEqual(get_hull_distance(), 0.10335952666666692)

    def test_get_competing_phases_for_BiTeCl(self):
        os.chdir(ROOT)
        os.chdir('BiTeCl')
        competing_phases = get_competing_phases()
        self.assertEqual(competing_phases, [(u'BiTeCl', u'mp-28944')])


if __name__ == '__main__':
    unittest.main()
