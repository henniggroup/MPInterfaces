from __future__ import print_function, division, unicode_literals

import os

from twod_materials.utils import *
from twod_materials.stability.startup import relax

from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.core.structure import Structure

from monty.serialization import loadfn

import numpy as np

import math

import twod_materials


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')
PACKAGE_PATH = '/'.join(PACKAGE_PATH.split('/')[:-2])

try:
    config_vars = loadfn(os.path.join(os.path.expanduser('~'), 'config.yaml'))
except:
    print('WARNING: No config.yaml file was found. please configure the '\
    'config.yaml and put it in your home directory.')
    # Still set them for testing purposes.
    config_vars = loadfn(os.path.join(PACKAGE_PATH, 'config.yaml'))
    VASP = config_vars['normal_binary']
if 'queue_system' in config_vars:
    QUEUE = config_vars['queue_system'].lower()
elif '/ufrc/' in os.getcwd():
    QUEUE = 'slurm'
elif '/scratch/' in os.getcwd():
    QUEUE = 'pbs'


def run_pbe_calculation(dim=2, submit=True, force_overwrite=False):
    """
    Setup and submit a normal PBE calculation for band structure along
    high symmetry k-paths.

    Args:
        dim (int): 2 for relaxing a 2D material, 3 for a 3D material.
        submit (bool): Whether or not to submit the job.
        force_overwrite (bool): Whether or not to overwrite files
            if an already converged vasprun.xml exists in the
            directory.
    """

    PBE_INCAR_DICT = {'EDIFF': 1e-6, 'IBRION': 2, 'ISIF': 3,
                      'ISMEAR': 1, 'NSW': 0, 'LVTOT': True, 'LVHAR': True,
                      'LORBIT': 1, 'LREAL': 'Auto', 'NPAR': 4,
                      'PREC': 'Accurate', 'LWAVE': True, 'SIGMA': 0.1,
                      'ENCUT': 500, 'ISPIN': 2}

    directory = os.getcwd().split('/')[-1]

    if not os.path.isdir('pbe_bands'):
        os.mkdir('pbe_bands')
    if force_overwrite or not is_converged('pbe_bands'):
        os.system('cp CONTCAR pbe_bands/POSCAR')
        if os.path.isfile('POTCAR'):
            os.system('cp POTCAR pbe_bands/')
        PBE_INCAR_DICT.update({'MAGMOM': get_magmom_string()})
        Incar.from_dict(PBE_INCAR_DICT).write_file('pbe_bands/INCAR')
        structure = Structure.from_file('POSCAR')
        kpath = HighSymmKpath(structure)
        Kpoints.automatic_linemode(20, kpath).write_file('pbe_bands/KPOINTS')
        os.chdir('pbe_bands')
        if dim == 2:
            remove_z_kpoints()
        if QUEUE == 'pbs':
            write_pbs_runjob(directory, 1, 16, '800mb', '6:00:00', VASP)
            submission_command = 'qsub runjob'

        elif QUEUE == 'slurm':
            write_slurm_runjob(directory, 16, '800mb', '6:00:00', VASP)
            submission_command = 'sbatch runjob'

        if submit:
            os.system(submission_command)

        os.chdir('../')


def run_hse_prep_calculation(dim=2, submit=True):
    """
    Submits a quick static calculation to calculate the IBZKPT
    file using a smaller number of k-points (200/atom instead of
    1000/atom). The other outputs from this calculation are
    essentially useless.

    Args:
        dim (int): 2 for relaxing a 2D material, 3 for a 3D material.
        submit (bool): Whether or not to submit the job.
    """

    if not os.path.isdir('hse_prep'):
        os.mkdir('hse_prep')
    os.chdir('hse_prep')
    os.system('cp ../CONTCAR ./POSCAR')
    if os.path.isfile('../POTCAR'):
        os.system('cp POTCAR .')
    relax(dim=2, submit=False)
    incar_dict = Incar.from_file('INCAR').as_dict()
    incar_dict.update({'NSW': 0, 'NELM': 1, 'LWAVE': False, 'LCHARG': False,
                       'LAECHG': False})
    Incar.from_dict(incar_dict).write_file('INCAR')

    Kpoints.automatic_density(
        Structure.from_file('POSCAR'), 200
    ).write_file('KPOINTS')

    if dim == 2:
        kpts_lines = open('KPOINTS').readlines()

        with open('KPOINTS', 'w') as kpts:
            for line in kpts_lines[:3]:
                kpts.write(line)
            kpts.write(kpts_lines[3].split()[0] + ' '
                       + kpts_lines[3].split()[1] + ' 1')

    if QUEUE == 'pbs':
        write_pbs_runjob('{}_prep'.format(
            os.getcwd().split('/')[-2]), 1, 16, '800mb', '6:00:00', VASP)
        submission_command = 'qsub runjob'

    elif QUEUE == 'slurm':
        write_slurm_runjob('{}_prep'.format(
            os.getcwd().split('/')[-2]), 16, '800mb', '6:00:00', VASP)
        submission_command = 'sbatch runjob'

    if submit:
        os.system(submission_command)

    os.chdir('../')


def run_hse_calculation(dim=2, submit=True, force_overwrite=False,
                        destroy_prep_directory=False):
    """
    Setup/submit an HSE06 calculation to get an accurate band structure.
    Requires a previous IBZKPT from a standard DFT run. See
    http://cms.mpi.univie.ac.at/wiki/index.php/Si_bandstructure for more
    details.

    Args:
        dim (int): 2 for relaxing a 2D material, 3 for a 3D material.
        submit (bool): Whether or not to submit the job.
        force_overwrite (bool): Whether or not to overwrite files
            if an already converged vasprun.xml exists in the
            directory.
        destroy_prep_directory (bool): whether or not to remove
            (rm -r) the hse_prep directory, if it exists. This
            can help you to automatically clean up and save space.
    """

    HSE_INCAR_DICT = {'LHFCALC': True, 'HFSCREEN': 0.2, 'AEXX': 0.25,
                      'ALGO': 'D', 'TIME': 0.4, 'NSW': 0,
                      'LVTOT': True, 'LVHAR': True, 'LORBIT': 11,
                      'LWAVE': True, 'NPAR': 8, 'PREC': 'Accurate',
                      'EDIFF': 1e-4, 'ENCUT': 450, 'ICHARG': 2, 'ISMEAR': 1,
                      'SIGMA': 0.1, 'IBRION': 2, 'ISIF': 3, 'ISPIN': 2}

    if not os.path.isdir('hse_bands'):
        os.mkdir('hse_bands')
    if force_overwrite or not is_converged('hse_bands'):
        os.chdir('hse_bands')
        os.system('cp ../CONTCAR ./POSCAR')
        if os.path.isfile('../POTCAR'):
            os.system('cp ../POTCAR .')
        HSE_INCAR_DICT.update({'MAGMOM': get_magmom_string()})
        Incar.from_dict(HSE_INCAR_DICT).write_file('INCAR')

        # Re-use the irreducible brillouin zone KPOINTS from a
        # previous standard DFT run.
        if os.path.isdir('../hse_prep'):
            ibz_lines = open('../hse_prep/IBZKPT').readlines()
            if destroy_prep_directory:
                os.system('rm -r ../hse_prep')
        else:
            ibz_lines = open('../IBZKPT').readlines()

        n_ibz_kpts = int(ibz_lines[1].split()[0])
        kpath = HighSymmKpath(Structure.from_file('POSCAR'))
        Kpoints.automatic_linemode(20, kpath).write_file('KPOINTS')
        if dim == 2:
            remove_z_kpoints()
        linemode_lines = open('KPOINTS').readlines()

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

        if QUEUE == 'pbs':
            write_pbs_runjob('{}_hsebands'.format(
                os.getcwd().split('/')[-2]), 2, 64, '1800mb', '50:00:00', VASP)
            submission_command = 'qsub runjob'

        elif QUEUE == 'slurm':
            write_slurm_runjob('{}_hsebands'.format(
                os.getcwd().split('/')[-2]), 64, '1800mb', '50:00:00', VASP)
            submission_command = 'sbatch runjob'

        if submit:
            os.system(submission_command)

        os.chdir('../')
