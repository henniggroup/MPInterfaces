import unittest
import os
from collections import defaultdict

from mpinterfaces.utils import *

from pymatgen import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


ROOT = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "mat2d", "stability", "tests")
)

class UtilsTest(unittest.TestCase):

    def test_is_converged(self):
        false_control = is_converged(ROOT)
        true_control = is_converged(os.path.join(ROOT, 'BiTeCl'))
        self.assertTrue(true_control)
        self.assertFalse(false_control)

    def test_ensure_vacuum_for_SiP(self):
        # Compound test for add_vacuum and get_spacing.
        os.chdir(ROOT)
        structure = Structure.from_file('POSCAR_SiP')
        structure = ensure_vacuum(structure, vacuum=15)
        self.assertAlmostEqual(get_spacing(structure), 15.0)

    def test_get_magmom_string_for_FeCl2(self):
        os.chdir(ROOT)
        structure = Structure.from_file('POSCAR_FeCl2')
        test_string = get_magmom_string(structure)
        self.assertEqual(test_string, u'1*6.0 2*0.5')

    def test_get_rotation_matrix(self):
        test_matrix = get_rotation_matrix([0, 0, 1], 2*np.pi)
        control_matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(test_matrix[i][j], control_matrix[i][j])

    def test_align_c_axis_for_non_aligned_structure(self):
        os.chdir(ROOT)
        structure = Structure.from_file('POSCAR_SF6')
        structure = align_axis(structure, 'c', (0, 0, 1))
        control_axis = [9.04099732e-13, -2.42627092e-13, 8.3829076073038635]
        for i in range(3):
            self.assertAlmostEqual(structure.lattice.matrix[2][i],
                                   control_axis[i])

    def test_align_c_axis_for_already_aligned_structure(self):
        os.chdir(ROOT)
        control_axis = [0, 0, 23.4186286267]
        structure = Structure.from_file('BiTeCl/POSCAR')
        structure = align_axis(structure, 'c', (0, 0, 1))
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
        test_file = open('KPOINTS')
        test_lines = test_file.readlines()
        control_file = open('circle_mesh_KPOINTS')
        control_lines = control_file.readlines()
        self.assertEqual(test_lines, control_lines)
        os.system('rm KPOINTS')
        test_file.close()
        control_file.close()

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
        test_file = open('KPOINTS')
        test_lines = test_file.readlines()
        print (test_lines)
        control_file = open('../BiTeCl_control/KPOINTS')
        control_lines = control_file.readlines()
        print (control_lines)
        self.assertEqual(test_lines, control_lines)
        os.system('rm KPOINTS')
        test_file.close()
        control_file.close()

    def test_get_run_cmmnd(self):
        os.chdir(os.path.join(ROOT, '../../../../'))
        QUEUE_SYSTEM='slurm'
        trial_output = get_run_cmmnd()
        correct_output = (defaultdict(None, {'account': None, 'mem': None, \
        'walltime': '10:00:00', 'nodes': 1, 'pre_rocket': None, 'job_name': None, \
        'ntasks': 16, 'email': None, 'rocket_launch': None}),None)
        self.assertEqual(trial_output, correct_output)



if __name__ == '__main__':
    unittest.main()
