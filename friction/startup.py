import os

import twod_materials.utils as utl

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar


def run_friction_calculations(directories, submit=True):

    for directory in directories:
        os.chdir(directory)
        utl.add_vacuum(3 - utl.get_spacing(), 0.8)
        structure = Structure.from_file('POSCAR')
        n_sites_per_layer = structure.num_sites
        structure.make_supercell([1, 1, 2])
        structure.to('POSCAR', 'POSCAR')
        utl.add_vacuum(12, 0.9)

        z_coords = []
        for site in structure.sites:
            z_coords.append(site.z)
        bottom_layer_max_height = sorted(z_coords)[n_sites_per_layer - 1]

        for x in range(10):
            for y in range(10):
                dir = '{}x{}'.format(x, y)

                # Copy input files
                os.system('cp INCAR {}/'.format(dir))
                os.system('cp KPOINTS {}/'.format(dir))
                os.system('cp POSCAR {}/'.format(dir))
                os.system('cp POTCAR {}/'.format(dir))
                os.system('cp vdw_kernel.bindat {}/'.format(dir))

                os.chdir(dir)
                incar_dict = Incar.from_file('INCAR').as_dict()
                incar_dict.update({'NSW': 0})
                Incar.from_dict(incar_dict).write_file('INCAR')

                # Shift the top layer
                poscar_lines = open('POSCAR').readlines()
                with open('POSCAR', 'w') as poscar:
                    for line in poscar_lines[:8]:
                        poscar.write(line)
                    for line in poscar_lines[8:8 + structure.num_sites]:
                        split_line = line.split()
                        if float(split_line[2]) > bottom_layer_max_height:
                            new_coords = [float(split_line[0]) + float(x)/10.0,
                                          float(split_line[1]) + float(y)/10.0,
                                          float(split_line[2])]
                            poscar.write(' '.join([str(i) for i in new_coords])
                                         + '\n')
                        else:
                            poscar.write(line)

                utl.write_runjob(dir, 1, 8, '400mb', '1:00:00', 'vasp')

                if submit:
                    os.system('qsub runjob')

                os.chdir('../')
        os.chdir('../')
