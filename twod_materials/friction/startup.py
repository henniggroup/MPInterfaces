import os

import math

import twod_materials.utils as utl

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar


def run_friction_calculations(spacing, submit=True):
    """
    Setup a 10x10 grid of static energy calculations to plot the Gamma
    surface between two layers of the 2D material.

    spacing: vertical distance (in angstroms) between the top of the
        bottom layer and the bottom of the top layer.
    """

    os.system('cp CONTCAR POSCAR')

    # Pad the bottom layer with 20 Angstroms of vacuum.
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
                           site.coords[2] + thickness + spacing]))

    for site in new_sites:
        structure.append(site[0], site[1], coords_are_cartesian=True)

    structure.to('POSCAR', 'POSCAR')

    n_divs_x = math.ceil(structure.lattice.a * 2.5)
    n_divs_y = math.ceil(structure.lattice.b * 2.5)

    for x in range(n_divs_x):
        for y in range(n_divs_y):
            dir = '{}x{}'.format(x, y)

            if not os.path.isdir(dir):
                os.mkdir(dir)

            # Copy input files
            os.system('cp INCAR {}/'.format(dir))
            os.system('cp KPOINTS {}/'.format(dir))
            os.system('cp POSCAR {}/'.format(dir))
            os.system('cp vdw_kernel.bindat {}/'.format(dir))

            os.chdir(dir)
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
                    new_coords = [float(split_line[0]) + float(x)/10.0,
                                  float(split_line[1]) + float(y)/10.0,
                                  float(split_line[2])]
                    poscar.write(' '.join([str(i) for i in new_coords])
                                 + '\n')

            utl.write_runjob(dir, 1, 8, '400mb', '1:00:00', 'vasp')

            if submit:
                os.system('qsub runjob')

            os.chdir('../')
    os.system('cp CONTCAR POSCAR')
