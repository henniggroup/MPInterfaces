import unittest

import os

from monty.serialization import loadfn

from pymatgen.core.structure import Structure
from pymatgen.matproj.rest import MPRester

import twod_materials
from twod_materials.electronic_structure.startup import run_hse_calculations
from twod_materials.electronic_structure.analysis import (get_band_structures,
                                                          plot_band_alignments)
from twod_materials.friction.startup import run_friction_calculations
from twod_materials.friction.analysis import plot_friction_surface
from twod_materials.pourbaix.startup import Calibrator
from twod_materials.pourbaix.analysis import Pourbaix2D
from twod_materials.stability.startup import relax, relax_competing_species
from twod_materials.stability.analysis import (get_competing_species,
                                               get_hull_distances)
from twod_materials.utils import (is_converged, add_vacuum, get_spacing,
                                  write_potcar, write_runjob)


PACKAGE_PATH = os.getcwd()
INCAR_DICT = {
    '@class': 'Incar', '@module': 'pymatgen.io.vasp.inputs', 'AGGAC': 0.0,
    'EDIFF': 1e-06, 'GGA': 'Bo', 'IBRION': 2, 'ISIF': 3, 'ISMEAR': 0,
    'LAECHG': True, 'LCHARG': True, 'LREAL': 'Auto', 'LUSE_VDW': True,
    'NPAR': 4, 'NSW': 50, 'PARAM1': 0.1833333333, 'PARAM2': 0.22,
    'PREC': 'High', 'SIGMA': 0.1
    }
KERNEL_PATH = os.path.join(PACKAGE_PATH, 'vdw_kernel.bindat')
MPR = MPRester(os.environ['MP_API'])
ION_DATA = loadfn(os.path.join(PACKAGE_PATH, 'pourbaix/ions.yaml'))
END_MEMBERS = loadfn(os.path.join(PACKAGE_PATH, 'pourbaix/end_members.yaml'))
ION_COLORS = loadfn(os.path.join(PACKAGE_PATH, 'pourbaix/ion_colors.yaml'))


class UtilsTest(unittest.TestCase):

    def test_is_converged_with_True_and_False_controls(self):
        with open('job.log', 'w') as joblog:
            joblog.write('reached required accuracy')
        true_control = is_converged('./')
        with open('job.log', 'w') as joblog:
            joblog.write(' ')
        false_control = is_converged('.')
        os.system('rm job.log')
        self.assertTrue(true_control)
        self.assertFalse(false_control)

    def test_add_vacuum_with_several_odd_structures(self):
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
