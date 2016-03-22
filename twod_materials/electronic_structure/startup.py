import os

from twod_materials.utils import (
    is_converged, write_pbs_runjob,
    write_slurm_runjob)

from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.core.structure import Structure


if '/ufrc/' in os.getcwd():
    HIPERGATOR = 2
elif '/scratch/' in os.getcwd():
    HIPERGATOR = 1

VASP = loadfn(os.path.join(os.path.expanduser('~'),
                           'config.yaml'))['normal_binary']


def remove_z_kpoints(filename='KPOINTS'):
    """
    Strips all paths from a linemode KPOINTS that include a z-component,
    since these are not relevant for 2D materials.
    """

    kpoint_lines = open(filename).readlines()
    with open(filename, 'w') as kpts:
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


def run_linemode_calculation(submit=True, force_overwrite=False):
    """
    Setup and submit a normal PBE calculation for band structure along
    high symmetry k-paths.
    """

    PBE_INCAR_DICT = {'EDIFF': 1e-6, 'IBRION': 2, 'ISIF': 3, 'ISMEAR': 1,
                      'NSW': 0, 'LVTOT': True, 'LVHAR': True, 'LORBIT': 11,
                      'LREAL': 'Auto', 'NPAR': 4, 'PREC': 'Accurate'
                      'LWAVE': True, 'SIGMA': 0.1, 'ENCUT': 500, 'ICHARG': 11}

    directory = os.getcwd().split('/')[-1]

    if not os.path.isdir('pbe_bands'):
        os.mkdir('pbe_bands')
    if force_overwrite or not is_converged('pbe_bands'):
        os.chdir('pbe_bands')
        os.system('cp ../CONTCAR ./POSCAR')
        os.system('cp ../POTCAR ./')
        os.system('cp ../CHGCAR ./')
        Incar.from_dict(PBE_INCAR_DICT).write_file('INCAR')
        structure = Structure.from_file('POSCAR')
        kpath = HighSymmKpath(structure)
        Kpoints.automatic_linemode(20, kpath).write_file('KPOINTS')
        remove_z_kpoints()
        if HIPERGATOR == 1:
            write_pbs_runjob(directory, 1, 16, '800mb', '6:00:00', VASP)
            submission_command = 'qsub runjob'

        elif HIPERGATOR == 2:
            write_slurm_runjob(directory, 16, '800mb', '6:00:00', VASP)
            submission_command = 'sbatch runjob'

        if submit:
            os.system(submission_command)

        os.chdir('../')


def run_hse_calculation(submit=True, force_overwrite=False):
    """
    Setup/submit an HSE06 calculation to get an accurate band structure.
    Requires a previous WAVECAR and IBZKPT from a standard DFT run. See
    http://cms.mpi.univie.ac.at/wiki/index.php/Si_bandstructure for more
    details.
    """

    HSE_INCAR_DICT = {'LHFCALC': True, 'HFSCREEN': 0.2, 'AEXX': 0.25,
                      'ALGO': 'D', 'TIME': 0.4, 'NSW': 0,
                      'LVTOT': True, 'LVHAR': True, 'LORBIT': 11,
                      'LWAVE': True, 'NPAR': 8, 'PREC': 'Accurate',
                      'EDIFF': 1e-6, 'ENCUT': 500, 'ISMEAR': 1, 'SIGMA': 0.1,
                      'IBRION': 2, 'ISIF': 3}

    if not os.path.isdir('hse_bands'):
        os.mkdir('hse_bands')
    if force_overwrite or not is_converged('hse_bands'):
        os.chdir('hse_bands')
        os.system('cp ../CONTCAR ./POSCAR')
        os.system('cp ../POTCAR ./POTCAR')
        os.system('cp ../WAVECAR ./')
        Incar.from_dict(HSE_INCAR_DICT).write_file('INCAR')

        # Re-use the irreducible brillouin zone KPOINTS from a
        # previous standard DFT run.
        ibz_lines = open('../IBZKPT').readlines()
        n_ibz_kpts = int(ibz_lines[1].split()[0])
        kpath = HighSymmKpath(Structure.from_file('POSCAR'))
        Kpoints.automatic_linemode(20, kpath).write_file('linemode_KPOINTS')
        remove_z_kpoints(filename='linemode_KPOINTS')
        linemode_lines = open('linemode_KPOINTS').readlines()

        abs_path = []
        i = 4
        while i < len(linemode_lines):
            start_kpt = linemode_lines[i].split()
            end_kpt = linemode_lines[i+1].split()
            increments = [
                (float(end_kpt[0]) - float(start_kpt[0])) / 20,
                (float(end_kpt[1]) - float(start_kpt[1])) / 20,
                (float(end_kpt[2]) - float(start_kpt[2])) / 20
            ]

            abs_path.append(start_kpt[:3] + ['0', start_kpt[4]])
            for n in range(1, 20):
                abs_path.append(
                    [str(float(start_kpt[0]) + increments[0] * n),
                     str(float(start_kpt[1]) + increments[1] * n),
                     str(float(start_kpt[2]) + increments[2] * n), '0']
                    )
            abs_path.append(end_kpt[:3] + ['0', end_kpt[4]])
            i += 3

        n_linemode_kpts = len(abs_path)

        with open('KPOINTS', 'w') as kpts:
            kpts.write('Automatically generated mesh\n')
            kpts.write('{}\n'.format(n_ibz_kpts + n_linemode_kpts))
            kpts.write('Reciprocal Lattice\n')
            for line in ibz_lines[3:]:
                kpts.write(line)
            for point in abs_path:
                kpts.write('{}\n'.format(' '.join(point)))

        if HIPERGATOR == 1:
            write_pbs_runjob('{}_hsebands'.format(
                os.getcwd().split('/')[-2]), 2, 64, '1800mb', '50:00:00', VASP)
            submission_command = 'qsub runjob'

        elif HIPERGATOR == 2:
            write_slurm_runjob('{}_hsebands'.format(
                os.getcwd().split('/')[-2]), 64, '1800mb', '50:00:00', VASP)
            submission_command = 'sbatch runjob'

        if submit:
            os.system(submission_command)

        os.system('rm linemode_KPOINTS')

        os.chdir('../')
