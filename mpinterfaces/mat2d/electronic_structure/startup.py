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
from mpinterfaces.utils import write_pbs_runjob, \
    write_slurm_runjob, is_converged, get_magmom_string, remove_z_kpoints

__author__ = "Michael Ashton, Joshua J. Gabriel"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton, Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

# TODO: document functions args and returns properly
# TODO: the run_* functions must be refactored to reduce code duplication

def get_markovian_path(points):
    """
     Calculates the shortest path connecting an array of 2D
     points. Returns the points in order on that path.

     Args:
         points

     Returns:
    """
    dist = lambda x, y: math.hypot(y[0] - x[0], y[1] - x[1])
    paths = [p for p in itertools.permutations(points)]
    path_distances = [sum(map(lambda x: dist(x[0], x[1]), zip(p[:-1], p[1:]))) for p in paths]
    min_index = np.argmin(path_distances)
    return paths[min_index]


def remove_z_kpoints_linemode(output='KPOINTS'):
    """
    Strips all k-points linemode KPOINTS that include a
    z-component, since these are not relevant for 2D materials.
    Then re-computes the markovian path between the remaining
    2D points and writes it over the KPOINTS file.

    Args:
        output (str)
    """

    kpoint_lines = open('KPOINTS').readlines()

    twod_kpoints = []
    labels = {}

    for i in range(4, len(kpoint_lines), 3):
         kpt_1 = kpoint_lines[i].split()
         kpt_2 = kpoint_lines[i+1].split()

         if float(kpt_1[2]) == 0.0 and [float(kpt_1[0]), float(kpt_1[1])] not in twod_kpoints:
             twod_kpoints.append([float(kpt_1[0]), float(kpt_1[1])])
             labels[kpt_1[4]] = [float(kpt_1[0]), float(kpt_1[1])]

         if float(kpt_2[2]) == 0.0 and [float(kpt_2[0]), float(kpt_2[1])] not in twod_kpoints:
             twod_kpoints.append([float(kpt_2[0]), float(kpt_2[1])])
             labels[kpt_2[4]] = [float(kpt_2[0]), float(kpt_2[1])]

    kpath = get_markovian_path(twod_kpoints)

    with open(output, 'w') as kpts:
         for line in kpoint_lines[:4]:
             kpts.write(line)

         for i in range(len(kpath)):
             label_1 = [l for l in labels if labels[l] == kpath[i]][0]
             if i == len(kpath) - 1:
                 kpt_2 = kpath[0]
                 label_2 = [l for l in labels if labels[l] == kpath[0]][0]
             else:
                 kpt_2 = kpath[i+1]
                 label_2 = [l for l in labels if labels[l] == kpath[i+1]][0]

             kpts.write(' '.join([str(kpath[i][0]), str(kpath[i][1]), '0.0 !', label_1]))
             kpts.write('\n')
             kpts.write(' '.join([str(kpt_2[0]), str(kpt_2[1]), '0.0 !', label_2]))
             kpts.write('\n\n')


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

    directory = os.path.basename(os.getcwd())

    if not os.path.isdir('pbe_bands'):
        os.mkdir('pbe_bands')

    if force_overwrite or not is_converged('pbe_bands'):
        shutil.copy("CONTCAR", "pbe_bands/POSCAR")
        if os.path.isfile('POTCAR'):
            shutil.copy("POTCAR", "pbe_bands")
        PBE_INCAR_DICT.update(
            {'MAGMOM': get_magmom_string(Structure.from_file('POSCAR'))})
        Incar.from_dict(PBE_INCAR_DICT).write_file('pbe_bands/INCAR')
        structure = Structure.from_file('POSCAR')
        kpath = HighSymmKpath(structure)
        Kpoints.automatic_linemode(20, kpath).write_file('pbe_bands/KPOINTS')
        os.chdir('pbe_bands')
        if dim == 2:
            remove_z_kpoints()
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
        shutil.copy('POTCAR', '.')
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
        HSE_INCAR_DICT.update(
            {'MAGMOM': get_magmom_string(Structure.from_file('POSCAR'))}
        )
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


def get_2D_hse_kpoints(struct_for_path, ibzkpth):
    """
    Args:
        struct_for_path: Structure from which linemode k-points will be generated.
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
