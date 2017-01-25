from __future__ import print_function, division, unicode_literals

import os

import math

import numpy as np

from monty.serialization import loadfn

import twod_materials.utils as utl

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar

import twod_materials
from twod_materials.stability.startup import get_magmom_string


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')
KERNEL_PATH = os.path.join(PACKAGE_PATH, 'vdw_kernel.bindat')

if '/ufrc/' in os.getcwd():
    HIPERGATOR = 2
elif '/scratch/' in os.getcwd():
    HIPERGATOR = 1

VASP = loadfn(os.path.join(os.path.expanduser('~'),
                           'config.yaml'))['normal_binary']


def run_gamma_calculations(submit=True, step_size=0.5):
    """
    Setup a 2D grid of static energy calculations to plot the Gamma
    surface between two layers of the 2D material.

    Step size is the distance between grid points in Angstroms.
    """

    if not os.path.isdir('friction'):
        os.mkdir('friction')
    os.chdir('friction')

    if not os.path.isdir('lateral'):
        os.mkdir('lateral')
    os.chdir('lateral')

    os.system('cp ../../CONTCAR POSCAR')

    # Pad the bottom layer with 20 Angstroms of vacuum.
    utl.add_vacuum(20 - utl.get_spacing(), 0.8)
    structure = Structure.from_file('POSCAR')
    n_sites_per_layer = structure.num_sites

    n_divs_x = int(math.ceil(structure.lattice.a / step_size))
    n_divs_y = int(math.ceil(structure.lattice.b / step_size))

    # Get the thickness of the material.
    max_height = max([site.coords[2] for site in structure.sites])
    min_height = min([site.coords[2] for site in structure.sites])
    thickness = max_height - min_height

    # Make a new layer.
    new_sites = []
    for site in structure.sites:
        new_sites.append((site.specie,
                          [site.coords[0], site.coords[1],
                           site.coords[2] + thickness + 3.5]))

    for site in new_sites:
        structure.append(site[0], site[1], coords_are_cartesian=True)

    #structure.get_sorted_structure().to('POSCAR', 'POSCAR')
    structure.to('POSCAR', 'POSCAR')

    for x in range(n_divs_x):
        for y in range(n_divs_y):
            dir = '{}x{}'.format(x, y)

            if not os.path.isdir(dir):
                os.mkdir(dir)

            # Copy input files
            os.chdir(dir)
            os.system('cp ../../../INCAR .')
            os.system('cp ../../../KPOINTS .')
            os.system('cp ../POSCAR .')
            os.system('cp {} .'.format(KERNEL_PATH))

            utl.write_potcar()
            incar_dict = Incar.from_file('INCAR').as_dict()
            incar_dict.update({'NSW': 0, 'LAECHG': False, 'LCHARG': False,
                               'LWAVE': False, 'MAGMOM': get_magmom_string()})
            incar_dict.pop('NPAR', None)
            Incar.from_dict(incar_dict).write_file('INCAR')

            # Shift the top layer
            poscar_lines = open('POSCAR').readlines()
            with open('POSCAR', 'w') as poscar:
                for line in poscar_lines[:8 + n_sites_per_layer]:
                    poscar.write(line)
                for line in poscar_lines[8 + n_sites_per_layer:]:
                    split_line = line.split()
                    new_coords = [
                        float(split_line[0]) + float(x)/float(n_divs_x),
                        float(split_line[1]) + float(y)/float(n_divs_y),
                        float(split_line[2])]
                    poscar.write(' '.join([str(i) for i in new_coords])
                                 + '\n')

            if HIPERGATOR == 1:
                utl.write_pbs_runjob(dir, 1, 4, '400mb', '1:00:00', VASP)
                submission_command = 'qsub runjob'

            elif HIPERGATOR == 2:
                utl.write_slurm_runjob(dir, 4, '400mb', '1:00:00', VASP)
                submission_command = 'sbatch runjob'

            if submit:
                os.system(submission_command)

            os.chdir('../')

    os.chdir('../../')


def run_normal_force_calculations(basin_and_saddle_dirs,
                                  spacings=np.arange(1.5, 4.25, 0.25),
                                  submit=True):
    """
    Set up and run static calculations of the basin directory
    and saddle directory (specified as a tuple) at specified
    interlayer spacings (by default, between 1.5 and 4 Angstroms)
    to get f_N and f_F.
    ex.
        run_normal_force_calculations(('0x0', '3x6'))
    or
        run_normal_force_calculations(get_basin_and_peak_locations())
    """

    spacings = [str(spc) for spc in spacings]

    os.chdir('friction')
    if not os.path.isdir('normal'):
        os.mkdir('normal')
    os.chdir('normal')

    for spacing in spacings:
        if not os.path.isdir(spacing):
            os.mkdir(spacing)

        for subdirectory in basin_and_saddle_dirs:

            os.system('cp -r ../lateral/{} {}/'.format(subdirectory, spacing))

            os.chdir('{}/{}'.format(spacing, subdirectory))
            structure = Structure.from_file('POSCAR')
            n_sites = len(structure.sites)
            top_layer = structure.sites[n_sites / 2:]
            bottom_of_top_layer = min(
                [z_coord for z_coord in [site.coords[2] for site in top_layer]])

            remove_indices = range(n_sites / 2, n_sites)

            structure.remove_sites(remove_indices)
            max_height = max([site.coords[2] for site in structure.sites])

            for site in top_layer:
                structure.append(
                    site.specie,
                    [site.coords[0],
                     site.coords[1],
                     site.coords[2] - bottom_of_top_layer
                     + max_height + float(spacing)],
                     coords_are_cartesian=True
                    )

            structure.to('POSCAR', 'POSCAR')

            if HIPERGATOR == 1:
                utl.write_pbs_runjob('{}_{}'.format(subdirectory, spacing), 1,
                    4, '400mb', '1:00:00', VASP)
                submission_command = 'qsub runjob'

            elif HIPERGATOR == 2:
                utl.write_slurm_runjob('{}_{}'.format(subdirectory, spacing), 4,
                    '400mb', '1:00:00', VASP)
                submission_command = 'sbatch runjob'

            if submit:
                os.system(submission_command)

            os.chdir('../../')

    os.chdir('../../')
