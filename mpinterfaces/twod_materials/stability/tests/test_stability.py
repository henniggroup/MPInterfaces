import os
import unittest

from mpinterfaces import MPINT_CONFIG
from mpinterfaces.twod_materials.stability import *

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "tests"))


class StartupTest(unittest.TestCase):

    @unittest.skipIf(not MPINT_CONFIG, "MPINT_CONFIG not set")
    def test_relax_creates_files(self):
        os.chdir(ROOT)
        os.chdir('BiTeCl')
        relax(submit=False, force_overwrite=True)
        for f in ['KPOINTS', 'INCAR', 'runjob']:
            self.assertTrue(os.path.isfile(f))
            os.system('rm {}'.format(f))


#class AnalysisTest(unittest.TestCase):

#    def test_get_hull_distance_for_BiTeCl(self):
#        os.chdir(ROOT)
#        os.chdir('BiTeCl')
#        self.assertEqual(get_hull_distance(), 0.10335952666666692)

#    def test_get_competing_phases_for_BiTeCl(self):
#        os.chdir(ROOT)
#        os.chdir('BiTeCl')
#        competing_phases = get_competing_phases()
#        self.assertEqual(competing_phases, [(u'BiTeCl', u'mp-28944')])


if __name__ == '__main__':
    unittest.main()
