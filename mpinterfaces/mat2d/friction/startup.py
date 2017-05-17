from __future__ import print_function, division, unicode_literals

import math
import os
import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar

from mpinterfaces import VASP_STD_BIN, VDW_KERNEL, QUEUE_SYSTEM
import mpinterfaces.utils as utl

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


# TODO: the run_* functions in mat2d subpackages can be merged and simplified, a lot
# of code duplcations

def run_gamma_calculations(submit=True, step_size=0.5):
    """
    Setup a 2D grid of static energy calculations to plot the Gamma
    surface between two layers of the 2D material. These calculations
    are run and stored in subdirectories under 'friction/lateral'.

    Args:
        submit (bool): Whether or not to submit the jobs.
        step_size (float): the distance between grid points in
            Angstroms.
    """

    if not os.path.isdir('friction'):
        os.mkdir('friction')
    os.chdir('friction')

    if not os.path.isdir('lateral'):
        os.mkdir('lateral')
    os.chdir('lateral')

    os.system('cp ../../CONTCAR POSCAR')

    # Pad the bottom layer with 20 Angstroms of vacuum.
    utl.ensure_vacuum(Structure.from_file('POSCAR'), 20)
    structure = Structure.from_file('POSCAR')
    n_sites_per_layer = structure.num_sites

    n_divs_x = int(math.ceil(structure.lattice.a / step_size))
    n_divs_y = int(math.ceil(structure.lattice.b / step_size))

    # Get the thickness of the material.
    max_height = max([site.coords[2] for site in structure.sites])
    min_height = min([site.coords[2] for site in structure.sites])
    thickness = max_height - min_height

    # Make a new layer.
    species, coords = [], []
    for site in structure.sites:
        # Original site
        species.append(site.specie)
        coords.append(site.coords)
        # New layer site
        species.append(site.specie)
        coords.append([site.coords[0], site.coords[1],
                       site.coords[2] + thickness + 3.5])

    Structure(structure.lattice, species, coords,
              coords_are_cartesian=True).to('POSCAR', 'POSCAR')

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
            if VDW_KERNEL:
                os.system('cp {} .'.format(VDW_KERNEL))

            # Shift the top layer
            structure = Structure.from_file("POSCAR")
            all_z_coords = [s.coords[2] for s in structure.sites]
            top_layer = [s for s in structure.sites if s.coords[2] > np.mean(all_z_coords)]
            structure.remove_sites([i for i, s in enumerate(structure.sites) if s in top_layer])
            for site in top_layer:
                structure.append(
                    site.specie,
                    [site.coords[0]+float(x)/float(n_divs_x),
                     site.coords[1]+float(y)/float(n_divs_y),
                     site.coords[2]], coords_are_cartesian=True
                )

            structure = structure.get_sorted_structure()
            structure.to("POSCAR", "POSCAR")
            utl.write_potcar()
            incar_dict = Incar.from_file('INCAR').as_dict()
            incar_dict.update({'NSW': 0, 'LAECHG': False, 'LCHARG': False,
                               'LWAVE': False, 'LVTOT': False,
                               'MAGMOM': utl.get_magmom_string(structure)})
            incar_dict.pop('NPAR', None)
            Incar.from_dict(incar_dict).write_file('INCAR')

            if QUEUE_SYSTEM == 'pbs':
                utl.write_pbs_runjob(dir, 1, 8, '1000mb', '2:00:00', VASP_STD_BIN)
                submission_command = 'qsub runjob'

            elif QUEUE_SYSTEM == 'slurm':
                utl.write_slurm_runjob(dir, 8, '1000mb', '2:00:00', VASP_STD_BIN)
                submission_command = 'sbatch runjob'

            if submit:
                os.system(submission_command)

            os.chdir('../')

    os.chdir('../../')


def run_normal_force_calculations(basin_and_saddle_dirs,
                                  spacings=np.arange(1.5, 4.25, 0.25),
                                  submit=True):
    """
    Set up and run static calculations of the basin directory and
    saddle directory at specified interlayer spacings to get f_N and
    f_F.

    Args:
        basin_and_saddle_dirs (tuple): Can be obtained by the
            get_basin_and_peak_locations() function under
            friction.analysis. For example,

            run_normal_force_calculations(('0x0', '3x6'))

            or

            run_normal_force_calculations(get_basin_and_peak_locations())

            will both work.
        spacings (tuple): list of interlayer spacings (in Angstroms, as floats)
            at which to run the calculations.
        submit (bool): Whether or not to submit the jobs.
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
            all_z_coords = [s.coords[2] for s in structure.sites]
            top_layer = [s for s in structure.sites if s.coords[2] >
                np.mean(all_z_coords)]
            bottom_of_top_layer = min([site.coords[2] for site in top_layer])

            remove_indices = [i for i, s in enumerate(structure.sites) if s in
                top_layer]
            structure.remove_sites(remove_indices)

            top_of_bottom_layer = max(
                [site.coords[2] for site in structure.sites]
            )

            for site in top_layer:
                structure.append(
                    site.specie,
                    [site.coords[0],
                     site.coords[1],
                     site.coords[2] - bottom_of_top_layer
                     + top_of_bottom_layer + float(spacing)],
                    coords_are_cartesian=True)

            structure = structure.get_sorted_structure()
            structure.to('POSCAR', 'POSCAR')
            utl.write_potcar()
            incar_dict = Incar.from_file('INCAR').as_dict()
            incar_dict.update({"MAGMOM": utl.get_magmom_string(structure)})
            Incar.from_dict(incar_dict).write_file("INCAR")

            if QUEUE_SYSTEM == 'pbs':
                utl.write_pbs_runjob('{}_{}'.format(
                    subdirectory, spacing), 1, 8, '1000mb', '2:00:00',
                    VASP_STD_BIN)
                submission_command = 'qsub runjob'

            elif QUEUE_SYSTEM == 'slurm':
                utl.write_slurm_runjob('{}_{}'.format(
                    subdirectory, spacing), 8, '1000mb', '2:00:00',
                    VASP_STD_BIN)
                submission_command = 'sbatch runjob'

            if submit:
                os.system(submission_command)

            os.chdir('../../')

    os.chdir('../../')
