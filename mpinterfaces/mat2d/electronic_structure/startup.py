from __future__ import print_function, division, unicode_literals

import itertools
import math
import numpy as np
import os
import shutil
import subprocess

from pymatgen import Structure
from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.symmetry.bandstructure import HighSymmKpath

from mpinterfaces import VASP_STD_BIN, QUEUE_SYSTEM
from mpinterfaces.mat2d.stability import relax
from mpinterfaces.utils import write_pbs_runjob, get_markovian_path,\
    write_slurm_runjob, is_converged, get_magmom_string, remove_z_kpoints

__author__ = "Michael Ashton, Joshua J. Gabriel"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton, Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

# TODO: document functions args and returns properly
# TODO: the run_* functions must be refactored to reduce code duplication


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

    PBE_INCAR_DICT = {'EDIFF': 1e-6, 'IBRION': 2, 'ICHARG': 2, 'ISIF': 3,
                      'ISMEAR': 1, 'NSW': 0, 'LVTOT': True, 'LVHAR': True,
                      'LORBIT': 1, 'LREAL': 'Auto', 'NPAR': 4,
                      'PREC': 'Accurate', 'LWAVE': True, 'SIGMA': 0.1,
                      'ENCUT': 500, 'ISPIN': 2}

    directory = os.path.basename(os.getcwd())

    if not os.path.isdir('pbe_bands'):
        os.mkdir('pbe_bands')

    if force_overwrite or not is_converged('pbe_bands'):
        shutil.copy("CONTCAR", "pbe_bands/POSCAR")
        structure = Structure.from_file("pbe_bands/POSCAR")
        if os.path.isfile("POTCAR"):
          shutil.copy("POTCAR", "pbe_bands")
        PBE_INCAR_DICT.update(
            {'MAGMOM': get_magmom_string(structure)})
        Incar.from_dict(PBE_INCAR_DICT).write_file('pbe_bands/INCAR')

        os.chdir('pbe_bands')
        write_band_structure_kpoints(structure, dim=dim)

        if QUEUE_SYSTEM == 'pbs':
            write_pbs_runjob(directory, 1, 16, '800mb', '6:00:00', VASP_STD_BIN)
            submission_command = 'qsub runjob'

        elif QUEUE_SYSTEM == 'slurm':
            write_slurm_runjob(directory, 16, '800mb', '6:00:00', VASP_STD_BIN)
            submission_command = 'sbatch runjob'

        if submit:
            _ = subprocess.check_output(submission_command.split())

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
    shutil.copy('../CONTCAR',  'POSCAR')
    if os.path.isfile('../POTCAR'):
        shutil.copy('../POTCAR', '.')
    relax(dim=2, submit=False)
    incar_dict = Incar.from_file('INCAR').as_dict()
    incar_dict.update({'NSW': 0, 'NELM': 1, 'LWAVE': False, 'LCHARG': False,
                       'LAECHG': False})
    Incar.from_dict(incar_dict).write_file('INCAR')

    Kpoints.automatic_density(
        Structure.from_file('POSCAR'), 200).write_file('KPOINTS')

    if dim == 2:
        kpts_lines = open('KPOINTS').readlines()

        with open('KPOINTS', 'w') as kpts:
            for line in kpts_lines[:3]:
                kpts.write(line)
            kpts.write(kpts_lines[3].split()[0] + ' '
                       + kpts_lines[3].split()[1] + ' 1')

    if QUEUE_SYSTEM == 'pbs':
        write_pbs_runjob('{}_prep'.format(
            os.getcwd().split('/')[-2]), 1, 16, '800mb', '6:00:00', VASP_STD_BIN)
        submission_command = 'qsub runjob'

    elif QUEUE_SYSTEM == 'slurm':
        write_slurm_runjob('{}_prep'.format(
            os.getcwd().split('/')[-2]), 16, '800mb', '6:00:00', VASP_STD_BIN)
        submission_command = 'sbatch runjob'

    if submit:
        _ = subprocess.check_output(submission_command.split())

    os.chdir('../')


def run_hse_calculation(dim=2, submit=True, force_overwrite=False):
    """
    Setup/submit an HSE06 calculation to get an accurate band structure.
    See http://cms.mpi.univie.ac.at/wiki/index.php/Si_bandstructure for
    more details.

    Args:
        dim (int): 2 for relaxing a 2D material, 3 for a 3D material.
        submit (bool): Whether or not to submit the job.
        force_overwrite (bool): Whether or not to overwrite files
            if an already converged vasprun.xml exists in the
            directory.
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
        structure = Structure.from_file("POSCAR")
        if os.path.isfile('../POTCAR'):
            os.system('cp ../POTCAR .')
        HSE_INCAR_DICT.update(
            {'MAGMOM': get_magmom_string(structure)}
        )
        Incar.from_dict(HSE_INCAR_DICT).write_file('INCAR')

        # Re-use the irreducible brillouin zone KPOINTS from a
        # previous standard DFT run.
        write_band_structure_kpoints(structure, dim=dim)

        if QUEUE_SYSTEM == 'pbs':
            write_pbs_runjob('{}_hsebands'.format(
                os.getcwd().split('/')[-2]), 2, 64, '1800mb', '50:00:00', VASP_STD_BIN)
            submission_command = 'qsub runjob'

        elif QUEUE_SYSTEM == 'slurm':
            write_slurm_runjob('{}_hsebands'.format(
                os.getcwd().split('/')[-2]), 64, '1800mb', '50:00:00', VASP_STD_BIN)
            submission_command = 'sbatch runjob'

        if submit:
            _ = subprocess.check_output(submission_command.split())

        os.chdir('../')


def write_band_structure_kpoints(structure, n_kpts=20, dim=2,
                                 ibzkpt_path="../"):
    """
    Writes a KPOINTS file for band structure calculations. Does
    not use the typical linemode syntax for NSCF calculations,
    but uses the IBZKPT + high-symmetry path syntax described in
    http://cms.mpi.univie.ac.at/wiki/index.php/Si_bandstructure
    so that SCF calculations can be performed. This is more
    reliable than re-using the CHGCAR from a previous run, which
    often results in "dimensions on the CHGCAR are different"
    errors in VASP.

    Args:
        structure (Structure): structure for determining k-path
        n_kpts (int): number of divisions along high-symmetry lines
        dim (int): 2 for a 2D material, 3 for a 3D material.
        ibzkpt_path (str): location of IBZKPT file. Defaults to one
            directory up.
    """

    ibz_lines = open(os.path.join(ibzkpt_path, "IBZKPT")).readlines()
    for i, line in enumerate(ibz_lines):
        if "Tetrahedra" in line:
            ibz_lines = ibz_lines[:i]
            break

    n_ibz_kpts = int(ibz_lines[1].split()[0])
    kpath = HighSymmKpath(structure)
    Kpoints.automatic_linemode(n_kpts, kpath).write_file('KPOINTS')
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


def get_2D_hse_kpoints(struct_for_path, ibzkpth):
    """
    Args:
        struct_for_path: Structure from which linemode k-points will
            be generated.
        ibzkpth:

    Returns:
        the Kpoints file object in the form of a string
              ready for execution by MPInterfaces
              calibrate objects
    """
    # Read IBZKPT from prep step
    ibz_lines = open(ibzkpth).readlines()
    n_ibz_kpts = int(ibz_lines[1].split()[0])

    # Read linemode KPOINTs from the dict (makes sure it is Kpoints
    # file with only 20 per atom for the optimized settings
    # Kpoints.from_dict(kpoint_dict).write_file('linemode_KPOINTS')
    kpath = HighSymmKpath(struct_for_path)
    Kpoints.automatic_linemode(20, kpath).write_file('KPOINTS_linemode')
    remove_z_kpoints_linemode()
    linemode_lines = open('KPOINTS_linemode').readlines()

    # put them together
    abs_path = []
    for i in range(4, len(linemode_lines), 3):
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

    n_linemode_kpts = len(abs_path)

    # write out the kpoints file and return the object

    Kpoints_hse_file = '\n'.join(
        ['Automatically generated mesh',
         '{}'.format(n_ibz_kpts + n_linemode_kpts),
         'Reciprocal Lattice',
         '{}'.format(str(''.join([line for line in ibz_lines[3:]])))]) + \
                       '{}'.format(str('\n'.join(
                           [' '.join(point) for point in abs_path])))

    ## can be used for test print out
    # with open('KPOINTS_HSE', 'w') as kpts:
    #        kpts.write('Automatically generated mesh\n')
    #        kpts.write('{}\n'.format(n_ibz_kpts + n_linemode_kpts))
    #        kpts.write('Reciprocal Lattice\n')
    #        for line in ibz_lines[3:]:
    #            kpts.write(line)
    #        for point in abs_path:
    #            kpts.write('{}\n'.format(' '.join(point)))

    return Kpoints_hse_file


def get_2D_incar_hse_prep(incar_dict):
    """
    linker for prep calculation

    Args:
        incar_dict (dict)

    Returns:
        dict: incar dict
    """
    print('updating INCAR for prep calculation ')
    INCAR_PREP = {'NSW': 0,
                  'NELM': 1,
                  'LWAVE': False,
                  'LCHARG': False,
                  'LAECHG': False}
    incar_dict.update(INCAR_PREP)
    return incar_dict


def get_2D_incar_hse(incar_dict):
    """
    linker function to complete the HSE input deck to MPInterfaces

    Args:
        incar_dict (dict)

    Returns:
        dict: incar dict
    """
    HSE_INCAR_DICT = {'LHFCALC': True, 'HFSCREEN': 0.2, 'AEXX': 0.25,
                      'ALGO': 'D', 'TIME': 0.4, 'NSW': 0, 'NELM': 75,
                      'LVTOT': True, 'LVHAR': True, 'LORBIT': 11,
                      'LWAVE': False, 'NPAR': 8, 'PREC': 'Accurate',
                      'EDIFF': 1e-4, 'ENCUT': 450, 'ICHARG': 2, 'ISMEAR': 1,
                      'SIGMA': 0.1, 'IBRION': 2, 'ISIF': 3, 'ISPIN': 2}
    incar_dict.update(HSE_INCAR_DICT)
    return incar_dict
