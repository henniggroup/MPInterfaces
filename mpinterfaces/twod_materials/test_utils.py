import unittest

import twod_materials
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath

from mpinterfaces.twod_materials.utils import *

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__version__ = "1.6"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')
ROOT = os.path.join(PACKAGE_PATH, 'stability/tests')


class UtilsTest(unittest.TestCase):


    def test_is_converged(self):
        false_control = is_converged(ROOT)
        true_control = is_converged(os.path.join(ROOT, 'BiTeCl'))
        self.assertTrue(true_control)
        self.assertFalse(false_control)


    def test_add_vacuum_and_get_spacing_with_odd_structures(self):
        os.chdir(ROOT)
        structure = Structure.from_file('POSCAR_SiP')
        top_layer = []
        for i in range(len(structure.sites)):
            if structure.sites[i].c > 0.5:
                top_layer.append(i)
        structure.remove_sites(top_layer)

        structure.to('POSCAR', 'POSCAR')
        add_vacuum(15 - get_spacing(), 0.5)
        self.assertTrue(14.9 < get_spacing() < 15.1)
        os.system('rm POSCAR')


    def test_get_magmom_string_for_FeCl2(self):
        os.chdir(ROOT)
        os.system('cp POSCAR_FeCl2 POSCAR')
        self.assertEqual(get_magmom_string(), u'1*6.0 2*0.5')
        os.system('rm POSCAR')


    def test_get_rotation_matrix(self):
        test_matrix = get_rotation_matrix((0, 0, 1), 2*np.pi)
        control_matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(test_matrix[i][j], control_matrix[i][j])


    def test_align_c_axis_for_non_aligned_structure(self):
        os.chdir(ROOT)
        structure = Structure.from_file('POSCAR_SF6')
        structure = align_c_axis_along_001(structure)
        control_axis = [9.04099732e-13, -2.42627092e-13, 8.3829076073038635]
        for i in range(3):
            self.assertAlmostEqual(structure.lattice.matrix[2][i],
                                   control_axis[i])


    def test_align_c_axis_for_already_aligned_structure(self):
        os.chdir(ROOT)
        control_axis = [0, 0, 23.4186286267]
        structure = Structure.from_file('BiTeCl/POSCAR')
        structure = align_c_axis_along_001(structure)
        for i in range(3):
            self.assertTrue(abs(
                structure.lattice.matrix[2][i] - control_axis[i]
                ) < 0.0001)


    def test_get_structure_type_for_conventional_material(self):
        os.chdir(ROOT)
        structure = Structure.from_file('POSCAR_Fe')
        test_type = get_structure_type(structure)
        self.assertEqual(test_type, 'conventional')


    def test_get_structure_type_for_layered_material(self):
        os.chdir(ROOT)
        structure = Structure.from_file('POSCAR_FeCl2')
        test_type = get_structure_type(structure)
        self.assertEqual(test_type, 'layered')


    def test_get_structure_type_for_1D_material(self):
        pass  # I still need to find one of these...


    def test_get_structure_type_for_0D_material(self):
        os.chdir(ROOT)
        structure = Structure.from_file('POSCAR_O2')
        test_type = get_structure_type(structure)
        self.assertEqual(test_type, 'molecular')


    def test_write_circle_mesh_kpoints(self):
        os.chdir(ROOT)
        write_circle_mesh_kpoints()
        test_lines = open('KPOINTS').readlines()
        control_lines = open('circle_mesh_KPOINTS').readlines()
        self.assertEqual(test_lines, control_lines)
        os.system('rm KPOINTS')


    def test_get_markovian_path(self):
        points = ((0, 0), (1, 1), (1, 0), (0, 1))
        control_points = ((0, 0), (1, 0), (1, 1), (0, 1))
        test_points = get_markovian_path(points)
        for i in range(len(control_points)):
            self.assertEqual(test_points[i], control_points[i])


    def test_remove_z_kpoints(self):
        os.chdir(os.path.join(ROOT, 'BiTeCl'))
        structure = Structure.from_file('POSCAR')
        kpath = HighSymmKpath(structure)
        Kpoints.automatic_linemode(20, kpath).write_file('KPOINTS')
        remove_z_kpoints()
        test_lines = open('KPOINTS').readlines()
        control_lines = open('../BiTeCl_control/KPOINTS').readlines()
        self.assertEqual(test_lines, control_lines)
        os.system('rm KPOINTS')


if __name__ == '__main__':
    unittest.main()
