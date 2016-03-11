import os

import math

import twod_materials.utils as utl

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar


if '/ufrc/' in os.getcwd():
    HIPERGATOR = 2
elif '/scratch/' in os.getcwd():
    HIPERGATOR = 1


def run_friction_calculations(submit=True):
    """
    Setup a 3D grid of static energy calculations to plot the Gamma
    surface between two layers of the 2D material, and obtain values of
    f_N for deriving friction coefficients.
    """

    # Pad the bottom layer with 20 Angstroms of vacuum.
    utl.add_vacuum(20 - utl.get_spacing(), 0.8)
    structure = Structure.from_file('POSCAR')
    n_sites_per_layer = structure.num_sites

    n_divs_x = int(math.ceil(structure.lattice.a * 2.5))
    n_divs_y = int(math.ceil(structure.lattice.b * 2.5))

    if not os.path.isdir('friction'):
        os.mkdir('friction')
    os.chdir('friction')

    os.system('cp ../CONTCAR POSCAR')
    utl.add_vacuum(20 - utl.get_spacing(), 0.8)
    structure = Structure.from_file('POSCAR')
    n_sites_per_layer = structure.num_sites

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

    structure.to('POSCAR', 'POSCAR')

    for x in range(n_divs_x):
        for y in range(n_divs_y):
            dir = '{}x{}'.format(x, y)

            if not os.path.isdir(dir):
                os.mkdir(dir)

            # Copy input files
            os.chdir(dir)
            print os.getcwd()
            os.system('cp ../../INCAR .')
            os.system('cp ../../KPOINTS .')
            os.system('cp ../POSCAR .')
            os.system('cp ../../vdw_kernel.bindat .')

            utl.write_potcar()
            incar_dict = Incar.from_file('INCAR').as_dict()
            incar_dict.update({'NSW': 0, 'LAECHG': False, 'LCHARG': False,
                               'LWAVE': False})
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
                utl.write_pbs_runjob(dir, 1, 4, '400mb', '1:00:00', 'vasp')
                submission_command = 'qsub runjob'

            elif HIPERGATOR == 2:
                utl.write_slurm_runjob(dir, 4, '400mb', '1:00:00', 'vasp')
                submission_command = 'sbatch runjob'

            if submit:
                os.system(submission_command)

            os.chdir('../')

    os.chdir('../')
