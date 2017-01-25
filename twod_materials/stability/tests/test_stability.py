import unittest

import os

from twod_materials.stability.analysis import (get_competing_phases,
                                               get_hull_distances)


ROOT = os.getcwd()
os.chdir('twod_materials/stability/tests')


class AnalysisTest(unittest.TestCase):

    def test_get_hull_distances_for_BiTeCl(self):
        self.assertEqual(get_hull_distances(['BiTeCl']),
                         {u'BiTeCl': 0.10335952666666692})

    def test_get_competing_phases_for_SiP(self):
        competing_phases = get_competing_phases(['BiTeCl'])
        self.assertEqual(competing_phases, [(u'BiTeCl', u'mp-28944')])

if __name__ == '__main__':
    unittest.main()
    os.chdir(ROOT)
