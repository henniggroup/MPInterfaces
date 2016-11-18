from __future__ import print_function, division, unicode_literals

import os

from twod_materials.utils import (
    is_converged, write_pbs_runjob,
    write_slurm_runjob)
from twod_materials.stability.startup import get_magmom_string, relax

from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.core.structure import Structure

from monty.serialization import loadfn

import numpy as np

import itertools as it

import math


if '/ufrc/' in os.getcwd():
    HIPERGATOR = 2
elif '/scratch/' in os.getcwd():
    HIPERGATOR = 1
else:
    HIPERGATOR = None

VASP = loadfn(os.path.join(os.path.expanduser('~'),
                           'config.yaml'))['normal_binary']


def write_circle_mesh_kpoints(center=(0, 0, 0), dimension=0.1,
                              resolution=20):
    """
    Create a circular mesh of k-points centered around a specific
    k-point (defaults to Gamma). `dimension` and `resolution`
    specify how large and how fine the grid will be. Non-circular
    meshes are not supported, but shouldn't be too hard to code.
    All k-point weights are 1.
    """

    kpoints = []
    step = dimension / resolution

    for i in range(-resolution, resolution):
        for j in range(-resolution, resolution):
            if i**2 + j**2 <= resolution**2:
                kpoints.append([str(center[0]+step*i), str(center[1]+step*j),
                '0', '1'])
    with open('KPOINTS', 'w') as kpts:
        kpts.write('KPOINTS\n{}\ndirect\n'.format(len(kpoints)))
        for kpt in kpoints:
            kpts.write(' '.join(kpt))
            kpts.write('\n')


def get_markovian_path(points):
    """
    Calculates the shortest path connecting an array of 2D
    points. Returns the points in order on that path.
    """

    def dist(x,y):
        return math.hypot(y[0] - x[0], y[1] - x[1])

    paths = [p for p in it.permutations(points)]
    path_distances = [
        sum(map(lambda x: dist(x[0], x[1]), zip(p[:-1], p[1:]))) for p in paths
    ]
    min_index = np.argmin(path_distances)

    return paths[min_index]


def remove_z_kpoints(output='KPOINTS'):
    """
    Strips all k-points linemode KPOINTS that include a
    z-component, since these are not relevant for 2D materials.
    Then re-computes the markovian path between the remaining
    2D points and writes it over the KPOINTS file.
    """

    kpoint_lines = open('KPOINTS').readlines()

    twod_kpoints = []
    labels = {}
    i = 4

    while i < len(kpoint_lines):
        kpt_1 = kpoint_lines[i].split()
        kpt_2 = kpoint_lines[i+1].split()
        if float(kpt_1[2]) == 0.0 and [float(kpt_1[0]), float(kpt_1[1])] not in twod_kpoints:
            twod_kpoints.append(
                [float(kpt_1[0]), float(kpt_1[1])]
            )
            labels[kpt_1[4]] = [float(kpt_1[0]), float(kpt_1[1])]

        if float(kpt_2[2]) == 0.0 and [float(kpt_2[0]), float(kpt_2[1])] not in twod_kpoints:
            twod_kpoints.append(
                [float(kpt_2[0]), float(kpt_2[1])]
            )
            labels[kpt_2[4]] = [float(kpt_2[0]), float(kpt_2[1])]
        i += 3

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


def run_linemode_calculation(submit=True, force_overwrite=False, dim='2D'):
    """
    Setup and submit a normal PBE calculation for band structure along
    high symmetry k-paths.

    `dim`: Set to "3D" if you want a 3D band structure. The default
    behavior is to remove all k-points with a z-component.
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
        os.system('cp POTCAR pbe_bands/')
        PBE_INCAR_DICT.update({'MAGMOM': get_magmom_string()})
        Incar.from_dict(PBE_INCAR_DICT).write_file('pbe_bands/INCAR')
        structure = Structure.from_file('POSCAR')
        kpath = HighSymmKpath(structure)
        Kpoints.automatic_linemode(20, kpath).write_file('pbe_bands/KPOINTS')
        os.chdir('pbe_bands')
        if dim == '2D':
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


def run_hse_prep_calculation(submit=True, dim='2D'):
    """
    Submits a quick static calculation to calculate the IBZKPT
    file using a smaller number of k-points (200/atom instead of
    1000/atom). The other outputs from this calculation are
    essentially useless.

    `dim`: Set to "3D" to include 3D k-points. The default behavior
    is for 2D materials, and removes all k-points with z-components.
    """

    if not os.path.isdir('hse_prep'):
        os.mkdir('hse_prep')
    os.chdir('hse_prep')
    os.system('cp ../CONTCAR ./POSCAR')
    relax(submit=False)
    incar_dict = Incar.from_file('INCAR').as_dict()
    incar_dict.update({'NSW': 0, 'NELM': 1, 'LWAVE': False, 'LCHARG': False,
                       'LAECHG': False})
    Incar.from_dict(incar_dict).write_file('INCAR')

    Kpoints.automatic_density(
        Structure.from_file('POSCAR'), 200
    ).write_file('KPOINTS')

    if dim == '2D':
        kpts_lines = open('KPOINTS').readlines()

        with open('KPOINTS', 'w') as kpts:
            for line in kpts_lines[:3]:
                kpts.write(line)
            kpts.write(kpts_lines[3].split()[0] + ' '
                       + kpts_lines[3].split()[1] + ' 1')

    if submit and HIPERGATOR == 1:
        os.system('qsub runjob')

    elif submit and HIPERGATOR == 2:
        os.system('sbatch runjob')


def run_hse_calculation(submit=True, force_overwrite=False,
                        destroy_prep_directory=False, dim='2D'):
    """
    Setup/submit an HSE06 calculation to get an accurate band structure.
    Requires a previous WAVECAR and IBZKPT from a standard DFT run. See
    http://cms.mpi.univie.ac.at/wiki/index.php/Si_bandstructure for more
    details.

    destroy_prep_directory=True will `rm -r` the hse_prep directory, if
    it exists. This can help you to automatically clean up and save
    space.

    `dim`: set to "3D" for a 3D band structure. The default behavior,
    "2D", is to remove all k-points with a z-component.
    """

    HSE_INCAR_DICT = {'LHFCALC': True, 'HFSCREEN': 0.2, 'AEXX': 0.25,
                      'ALGO': 'D', 'TIME': 0.4, 'NSW': 0,
                      'LVTOT': True, 'LVHAR': True, 'LORBIT': 11,
                      'LWAVE': True, 'NPAR': 8, 'PREC': 'Accurate',
                      'EDIFF': 1e-6, 'ENCUT': 500, 'ICHARG': 2, 'ISMEAR': 1,
                      'SIGMA': 0.1, 'IBRION': 2, 'ISIF': 3, 'ISPIN': 2}

    if not os.path.isdir('hse_bands'):
        os.mkdir('hse_bands')
    if force_overwrite or not is_converged('hse_bands'):
        os.chdir('hse_bands')
        os.system('cp ../CONTCAR ./POSCAR')
        os.system('cp ../POTCAR ./')
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
        if dim == '2D':
            remove_z_kpoints(output='KPOINTS')
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

        os.chdir('../')
