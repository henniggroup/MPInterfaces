import os

from pymatgen.core.structure import Structure

from monty.serialization import loadfn


POTENTIAL_PATH = loadfn(
    os.path.join(os.path.expanduser('~'), 'config.yaml'))['potentials']


def is_converged(directory):

    try:
        if 'reached required accuracy' in open(
            '{}/job.log'.format(directory)
                ).read():
            return True
        else:
            return False
    except IOError:
        return False


def get_spacing(filename='POSCAR', cutoff=0.95):
    """
    Returns the interlayer spacing for a 2D material.
    """

    structure = Structure.from_file('POSCAR')

    lines = open(filename).readlines()
    c_axis = lines[4].split()
    lattice_parameter = lines[1].split()
    split_coords = [line.split() for line in lines[8:8+structure.num_sites]]
    z_coords = list()
    for coord in split_coords:
        z_coord = float(coord[2])
        if z_coord > cutoff:
            z_coord -= 1
        z_coords.append(z_coord)
    max_height = max([z_height for z_height in z_coords])
    min_height = min([z_height for z_height in z_coords])
    spacing = ((1.0 + min_height) - max_height) * float(c_axis[2])\
        * float(lattice_parameter[0])

    return spacing


def add_vacuum(delta, cut=0.9):
    '''
    Adds vacuum to a POSCAR.

    delta = vacuum thickness in Angstroms
    cut = height above which atoms will need to be fixed. Defaults to
    0.9.
    '''
    # Fix the POSCAR to put bottom atoms (if they are accidentally above
    # tolerance) at 0.0.

    structure = Structure.from_file('POSCAR')
    with open('POSCAR', 'r') as readcar:
        poscar_lines = readcar.readlines()
    n_atoms = structure.num_sites
    atom_lines = []
    for i in range(8, 8+n_atoms):
        atom_lines.append(poscar_lines[i].split())
    atom_line_2s = []
    for atom_line in atom_lines:
        atom_line_2s.append(float(atom_line[2]))
    fixable = False
    addables = []
    for atom_line_2 in atom_line_2s:
        if float(atom_line_2) > cut or\
                (float(atom_line_2) < 0.0 and 1.0 + float(atom_line_2) > cut):
            if float(atom_line_2) < 0.0 and 1.0 + float(atom_line_2) > cut:
                atom_line_2 = float(atom_line_2) + 1.0
            addables.append(atom_line_2)
            fixable = True
    if fixable:
        add_factor = 1.0 - min(addables)
    else:
        add_factor = 0.0
    new_atom_lines = []
    for atom_line in atom_lines:
        new_atom_line_2 = str(float(atom_line[2]) + add_factor)
        if float(new_atom_line_2) >= 1.0:
            new_atom_line_2 = str(float(new_atom_line_2) - 1.0)
        new_atom_lines.append('{} {} {}'.format(atom_line[0], atom_line[1],
                                                new_atom_line_2))
    with open('POSCAR', 'w') as writecar:
        for line in poscar_lines[0:8]:
            writecar.write(line)
        for new_atom_line in new_atom_lines:
            writecar.write('{}\n'.format(new_atom_line))

    # Open files and read in values from POSCAR
    old_poscar = open('POSCAR', 'r')
    new_poscar = open('new_POSCAR', 'w')
    oldlines = old_poscar.readlines()
    name = oldlines[0].split()[0]
    lattice_constant = oldlines[1].strip()
    a_lat_par = oldlines[2].split()
    b_lat_par = oldlines[3].split()
    c_lat_par = oldlines[4].split()
    elements = oldlines[5].split()
    stoichiometry = oldlines[6].split()
    coordinate_type = oldlines[7].strip()
    n_atoms = 0
    for item in stoichiometry:
        n_atoms += int(item)
    final_atom_line = n_atoms + 8

    # Elongate c-vector by delta

    save = float(c_lat_par[2])
    c_length = float(c_lat_par[2]) * float(lattice_constant)
    c_length_plus_delta = c_length + float(delta)
    c_lat_par[2] = c_length_plus_delta / float(lattice_constant)
    scalar = c_lat_par[2] / save

    # Create list of atom coordinates and adjust their z-coordinate on
    # the fly

    atoms = []
    for i in range(8, final_atom_line):
        atom = oldlines[i].split()
        atom[2] = float(atom[2]) / scalar
        atoms.append(atom)

    # Write updated values to new_POSCAR, copy it to old_POSCAR, then
    # close files and delete new_POSCAR

    new_poscar.write('{}\n'.format(name))
    new_poscar.write('{}\n'.format(lattice_constant))
    for item in a_lat_par:
        new_poscar.write('{} '.format(item))
    new_poscar.write('\n')
    for item in b_lat_par:
        new_poscar.write('{} '.format(item))
    new_poscar.write('\n')
    for item in c_lat_par:
        new_poscar.write('{} '.format(item))
    new_poscar.write('\n')
    for item in elements:
        new_poscar.write('{} '.format(item))
    new_poscar.write('\n')
    for item in stoichiometry:
        new_poscar.write('{} '.format(item))
    new_poscar.write('\n')
    new_poscar.write('{}\n'.format(coordinate_type))
    for item in atoms:
        new_poscar.write('{} {} {}\n'.format(item[0], item[1], item[2]))

    new_poscar.close()
    old_poscar.close()
    os.remove('POSCAR')

    new_poscar = open('new_POSCAR', 'r')
    new_lines = new_poscar.readlines()
    poscar = open('POSCAR', 'w')
    for line in new_lines:
        poscar.write(line)
    old_poscar.close()
    poscar.close()
    new_poscar.close()
    os.remove('new_POSCAR')


def write_potcar(types='None'):
    '''
    Writes a POTCAR file based on a list of types.

    types = list of same length as number of elements containing specifications
    for the kind of potential desired for each element. If no special potential
    is desired, just enter 'regular', or leave types = 'None'.
    (['pv', 'regular', '3'])
    '''
    poscar = open('POSCAR', 'r')
    lines = poscar.readlines()
    elements = lines[5].split()
    poscar.close()

    if types == 'None':
        types = []
        for i in range(len(elements)):
            types.append('regular')
    potentials = []
    for i in range(len(elements)):
        if types[i] == 'regular':
            pass
        else:
            elements[i] += '_{}'.format(types[i])

        # If specified pseudopotential doesn't exist, try other variations.
        if os.path.exists('{}/{}/POTCAR'.format(POTENTIAL_PATH, elements[i])):
            pass
        else:
            print 'Potential file for {} does not exist. Looking for best'\
                  'variation... '.format(elements[i])
            if types[i] == 'regular':
                length = 0
            else:
                length = len(types[i]) + 1
                elements[i] = elements[i][:-length]
            elements[i] += '_sv'
            if os.path.exists('{}/{}/POTCAR'.format(
                    POTENTIAL_PATH, elements[i])):
                print 'Found one! {} will work.'.format(elements[i])
            else:
                elements[i] = elements[i][:-3]
                elements[i] += '_pv'
                if os.path.exists('{}/{}/POTCAR'.format(
                        POTENTIAL_PATH, elements[i])):
                    print 'Found one! {} will work.'.format(elements[i])
                else:
                    elements[i] = elements[i][:-3]
                    elements[i] += '_3'
                    if os.path.exists('{}/{}/POTCAR'.format(
                            POTENTIAL_PATH, elements[i])):
                        print 'Found one! {} will work.'.format(elements[i])
                    else:
                        elements[i] = elements[i][:-2]
                        if os.path.exists('{}/{}/POTCAR'.format(
                                POTENTIAL_PATH, elements[i])):
                            print ('Found one! {} will '
                                   'work.'.format(elements[i]))
                        else:
                            print ('No pseudopotential found'
                                   ' for {}'.format(elements[i]))

    # Create paths, open files, and write files to POTCAR for each potential.
    for element in elements:
        potentials.append('{}/{}/POTCAR'.format(POTENTIAL_PATH, element))
    outfile = open('POTCAR', 'w')
    for potential in potentials:
        infile = open(potential)
        for line in infile:
            outfile.write(line)
        infile.close()
    outfile.close()


def write_runjob(name, nnodes, nprocessors, pmem, walltime, binary):
    '''
    writes a runjob based on a name, nnodes, nprocessors, walltime, and
    binary. Designed for runjobs on the Hennig group_list on HiperGator.
    '''
    runjob = open('runjob', 'w')
    runjob.write('#!/bin/sh\n')
    runjob.write('#PBS -N {}\n'.format(name))
    runjob.write('#PBS -o test.out\n')
    runjob.write('#PBS -e test.err\n')
    runjob.write('#PBS -r n\n')
    runjob.write('#PBS -l walltime={}\n'.format(walltime))
    runjob.write('#PBS -l nodes={}:ppn={}\n'.format(nnodes, nprocessors))
    runjob.write('#PBS -l pmem={}\n'.format(pmem))
    runjob.write('#PBS -W group_list=hennig\n\n')
    runjob.write('cd $PBS_O_WORKDIR\n\n')
    runjob.write('mpirun ~/bin/{} > job.log\n\n'.format(binary))
    runjob.write('echo \'Done.\'\n')
    runjob.close()
