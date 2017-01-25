import unittest

import os

from monty.serialization import loadfn

from pymatgen.matproj.rest import MPRester

from twod_materials.utils import is_converged, add_vacuum, get_spacing


PACKAGE_PATH = os.path.join(os.getcwd(), 'twod_materials')
INCAR_DICT = {
    '@class': 'Incar', '@module': 'pymatgen.io.vasp.inputs', 'AGGAC': 0.0,
    'EDIFF': 1e-06, 'GGA': 'Bo', 'IBRION': 2, 'ISIF': 3, 'ISMEAR': 1,
    'LAECHG': True, 'LCHARG': True, 'LREAL': 'Auto', 'LUSE_VDW': True,
    'NPAR': 4, 'NSW': 50, 'PARAM1': 0.1833333333, 'PARAM2': 0.22,
    'PREC': 'High', 'SIGMA': 0.1
    }
KERNEL_PATH = os.path.join(PACKAGE_PATH, 'vdw_kernel.bindat')
ION_DATA = loadfn(os.path.join(PACKAGE_PATH, 'pourbaix/ions.yaml'))
END_MEMBERS = loadfn(os.path.join(PACKAGE_PATH, 'pourbaix/end_members.yaml'))
ION_COLORS = loadfn(os.path.join(PACKAGE_PATH, 'pourbaix/ion_colors.yaml'))
try:
    MPR = MPRester(
        loadfn(os.path.join(os.path.expanduser('~'), 'config.yaml'))['mp_api']
        )
except IOError:
    try:
        MPR = MPRester(
            os.environ['MP_API']
            )
    except KeyError:
        raise ValueError('No Materials Project API key found. Please check'
                         ' that your ~/config.yaml contains the field'
                         ' mp_api: your_api_key')


class UtilsTest(unittest.TestCase):

    def test_is_converged_with_controls(self):
        false_control = is_converged('./')
        true_control = is_converged(os.path.join(PACKAGE_PATH,
                                                 'stability/tests/BiTeCl'))
        self.assertTrue(true_control)
        self.assertFalse(false_control)

    def test_add_vacuum_with_odd_structures(self):
        structure = MPR.get_structure_by_material_id('mp-2798')  # SiP
        top_layer = []
        for i in range(len(structure.sites)):
            if structure.sites[i].c > 0.5:
                top_layer.append(i)
        structure.remove_sites(top_layer)

        structure.to('POSCAR', 'POSCAR')
        add_vacuum(15 - get_spacing(), 0.5)
        self.assertTrue(14.9 < get_spacing() < 15.1)
        os.system('rm POSCAR')


if __name__ == '__main__':
    unittest.main()
