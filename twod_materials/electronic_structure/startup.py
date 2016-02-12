import os

from twod_materials.utils import is_converged, write_runjob

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.core.structure import Structure


HSE_INCAR_DICT = {}


def remove_z_kpoints():
    """
    Strips all paths from a linemode KPOINTS that include a z-component,
    since these are not relevant for 2D materials.
    """

    kpoint_lines = open('KPOINTS').readlines()
    with open('KPOINTS', 'w') as kpts:
        for line in kpoint_lines[:4]:
            kpts.write(line)
        i = 4
        while i < len(kpoint_lines):
            if (
                not float(kpoint_lines[i].split()[2]) and
                not float(kpoint_lines[i+1].split()[2])
                    ):
                kpts.write(kpoint_lines[i])
                kpts.write(kpoint_lines[i+1])
                kpts.write('\n')
            i += 3


def run_linemode_calculation(submit=True):
    """
    Setup and submit a normal PBE calculation for band structure along
    high symmetry k-paths.
    """

    directory = os.getcwd().split('/')[-1]
    if is_converged('.'):
        if not os.path.isdir('bandstructure'):
            os.mkdir('bandstructure')
        os.chdir('bandstructure')
        os.system('cp ../CONTCAR ./POSCAR')
        os.system('cp ../POTCAR ./')
        os.system('cp ../vdw_kernel.bindat ./')
        incar_dict = Incar.from_file('../INCAR').as_dict()
        incar_dict.update({'NSW': 0, 'ISMEAR': -5})
        Incar.from_dict(incar_dict).write_file('INCAR')
        structure = Structure.from_file('POSCAR')
        kpath = HighSymmKpath(structure)
        Kpoints.automatic_linemode(20, kpath).write_file('KPOINTS')
        remove_z_kpoints()
        write_runjob(directory, 1, 16, '600mb', '6:00:00', 'vasp')

        if submit:
            os.system('qsub runjob')

        os.chdir('../')


def run_hse_calculation(submit=True):
    """
    Setup/submit an HSE06 calculation to get an accurate band structure.
    Requires a previous WAVECAR and IBZKPT from a PBE run. See
    http://cms.mpi.univie.ac.at/wiki/index.php/Si_bandstructure for more
    details.
    """

    directory = os.getcwd().split('/')[-1]
    vasprun = Vasprun('vasprun.xml')
    band_gap = vasprun.get_band_structure().get_band_gap()
    kpath = []
    if band_gap['energy']:
        transition = band_gap['transition'].split('-')

        # Amount to increment along the kpath to result in 10
        # points total.
        increment = ((transition[1][0] - transition[0][0]) / 9,
                     (transition[1][1] - transition[0][1]) / 9,
                     (transition[1][2] - transition[0][2] / 9))

        for i in range(10):
            kpath.append((transition[0] + increment[0] * i),
                         (transition[1] + increment[1] * i),
                         (transition[2] + increment[2] * i))

        if not os.path.isdir('HSE'):
            os.mkdir('HSE')
        os.chdir('HSE')
        os.system('cp ../CONTCAR ./POSCAR')
        os.system('cp ../POTCAR ./POTCAR')
        os.system('cp ../vdw_kernel.bindat ./')
        Incar.from_dict(HSE_INCAR_DICT).write_file('INCAR')
        write_runjob(directory, 1, 32, '1200mb', '150:00:00', 'vasp')

        # Re-use the irreducible brillouin zone KPOINTS from a
        # previous GGA run.
        os.system('cp ../IBZKPT ./KPOINTS')
        kpoints_lines = open('KPOINTS').readlines()
        n_kpts = int(kpoints_lines[1].split()[0])
        with open('KPOINTS', 'w') as kpts:
            kpts.write(kpoints_lines[0])
            kpts.write('{}\n'.format(n_kpts + 10))
            for line in kpoints_lines[2:]:
                kpts.write(line)
            for kpoint in kpath:
                kpts.write('{}\n'.format(' '.join(kpoint)))

        if submit:
            os.system('qsub runjob')
