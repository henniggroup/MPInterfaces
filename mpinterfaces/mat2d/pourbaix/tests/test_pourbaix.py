import os
import unittest

from monty.serialization import loadfn

from mpinterfaces import MP_API
from mpinterfaces.mat2d.pourbaix import *
from mpinterfaces.mat2d.stability.startup import INCAR_DICT

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "tests"))


class StartupTest(unittest.TestCase):

    def test_get_experimental_formation_energies(self):
        test_energies = get_experimental_formation_energies()
        control_energies = {
            'Ba': -5.503182859746106,
            'Be': -6.086201714841334,
            'Bi': -4.746555698345675,
            'C': -3.895039017966623,
            'Ca': -6.346175831299141,
            'Cd': -2.3003950526623558
        }
        for elt in control_energies:
            self.assertEqual(test_energies[elt], control_energies[elt])

    @unittest.skipIf(MP_API == None, "API key not set")
#    def test_relax_references_for_Mo_and_S(self):
#        os.chdir(ROOT)
#        relax_references(['Mo_pv', 'S'], incar_dict=INCAR_DICT, submit=False)
#        for d in ['Mo', 'Mo/oxide', 'S', 'O']:
#            self.assertTrue(os.path.isdir(d))
#            for f in ['POSCAR', 'INCAR', 'KPOINTS', 'runjob']:
#                self.assertTrue(os.path.isfile('{}/{}'.format(d, f)))
#        for d in ['Mo', 'S', 'O']:
#            os.system('rm -r {}'.format(d))

    def test_get_corrections_for_Mo_Ta_and_W(self):
        os.chdir(os.path.join(ROOT, 'Mo_Ta_W_controls'))
        get_corrections(write_yaml=True)
        test_corrections = loadfn('ion_corrections.yaml')
        control_corrections = {'Mo': 0.379, 'Ta': 0.161, 'W': 0.635}
        test_mus = loadfn('chemical_potentials.yaml')
        control_mus = {'Mo': -8.885, 'Ta': -10.17, 'W': -11.179, 'O': -4.658}

        for elt in control_corrections:
            self.assertEqual(control_corrections[elt], test_corrections[elt])
        for elt in control_mus:
            self.assertEqual(control_mus[elt], test_mus[elt])
        os.system('rm ion_corrections.yaml')
        os.system('rm chemical_potentials.yaml')
        os.chdir(ROOT)


if __name__ == '__main__':
    unittest.main()
