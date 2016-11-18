from __future__ import print_function, division, unicode_literals

import os

from pymatgen.core.structure import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.outputs import Vasprun

from monty.serialization import loadfn

import numpy as np

import math

import twod_materials


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')

try:
    POTENTIAL_PATH = loadfn(
        os.path.join(os.path.expanduser('~'), 'config.yaml'))['potentials']
    USR = loadfn(os.path.join(os.path.expanduser('~'),
                 'config.yaml'))['username']

except IOError:
    try:
        POTENTIAL_PATH = os.environ['VASP_PSP_DIR']
        USR = os.environ['USERNAME']
    except KeyError:
        raise ValueError('No config.yaml file found. Please check'
                         ' that your config.yaml is located in ~/ and'
                         ' contains the field'
                         ' potentials: /path/to/your/POTCAR/files/')


def is_converged(directory):
    """
    Check if a relaxation has converged.
    """

    try:
        if Vasprun('{}/vasprun.xml'.format(directory)).converged:
            return True
        else:
            return False

    except:
        return False


def get_status(directory):
    """
    Return the state of job in a directory. Designed for use on
    HiperGator.

    'C': complete
    'R': running
    'Q': queued
    'E': error
    'H': hold
    'None': No job in this directory
    """

    os.system("qstat -f| grep -A 30 '{}' >> my_jobs.txt".format(USR))
    lines = open('my_jobs.txt').readlines()
    job_state = None
    for i in range(len(lines)):
        if 'Output_Path' in lines[i]:
            joined_line = ''.join([lines[i].strip(), lines[i+1].strip()])
            if directory in joined_line:
                for j in range(i, 0, -1):
                    if 'job_state' in lines[j]:
                        job_state = lines[j].split('=')[1].strip()
                        break
    os.system('rm my_jobs.txt')

    return job_state


def get_spacing(filename='POSCAR', cut=0.9):
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
        if z_coord > cut:
            z_coord -= 1
        z_coords.append(z_coord)
    max_height = max([z_height for z_height in z_coords])
    min_height = min([z_height for z_height in z_coords])
    spacing = ((1.0 + min_height) - max_height) * abs(float(c_axis[2]))\
        * float(lattice_parameter[0])

    return spacing


def get_rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation
    about the given axis by theta radians.
    Credit: http://stackoverflow.com/users/190597/unutbu
    """

    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def align_c_axis_along_001(structure):
    """
    Given a structure with a c-axis not along [001],
    returns the same structure rotated so that the c-axis is along
    the [001] direction. This is useful for vasp compiled with no
    z-axis relaxation.
    """

    c = structure.lattice._matrix[2]
    z = [0, 0, 1]
    axis = np.cross(c, z)
    if not(axis[0] == 0 and axis[1] == 0):
        theta = (np.arccos(np.dot(c, z) /
                 (np.linalg.norm(c) * np.linalg.norm(z))))
        R = get_rotation_matrix(axis, theta)
        rotation = SymmOp.from_rotation_and_translation(rotation_matrix=R)
        structure.apply_operation(rotation)
    return structure


def get_structure_type(structure, write_poscar_from_cluster=False):
    """
    Returns 'molecular', 'chain', 'layered', 'heterogeneous', or
    'conventional' to describe the periodicity of bonded clusters
    in a bulk structure.

    Args:
        write_poscar_from_cluster (bool): Set to True to write a
            POSCAR from the sites in the cluster.
    """

    # The conventional standard structure is much easier to work
    # with.

    structure = SpacegroupAnalyzer(
        structure).get_conventional_standard_structure()

    # Noble gases don't have well-defined bonding radii.
    if not len([e for e in structure.composition
            if e.symbol in ['He', 'Ne', 'Ar', 'Kr', 'Xe']]) == 0:
        type = 'noble gas'
    else:
        if len(structure.sites) < 45:
            structure.make_supercell(2)

        # Create a dict of sites as keys and lists of their
        # bonded neighbors as values.
        sites = structure.sites
        bonds = {}
        for site in sites:
            bonds[site] = []

        for i in range(len(sites)):
            site_1 = sites[i]
            for site_2 in sites[i+1:]:
                if (site_1.distance(site_2) <
                        float(Element(site_1.specie).atomic_radius
                        + Element(site_2.specie).atomic_radius) * 1.1):
                    bonds[site_1].append(site_2)
                    bonds[site_2].append(site_1)

        # Assimilate all bonded atoms in a cluster; terminate
        # when it stops growing.
        cluster_terminated = False
        while not cluster_terminated:
            original_cluster_size = len(bonds[sites[0]])
            for site in bonds[sites[0]]:
                bonds[sites[0]] += [
                    s for s in bonds[site] if s not in bonds[sites[0]]]
            if len(bonds[sites[0]]) == original_cluster_size:
                cluster_terminated = True

        original_cluster = bonds[sites[0]]

        if len(bonds[sites[0]]) == 0:  # i.e. the cluster is a single atom.
            type = 'molecular'
        elif len(bonds[sites[0]]) == len(sites): # i.e. all atoms are bonded.
            type = 'conventional'
        else:
            # If the cluster's composition is not equal to the
            # structure's overall composition, it is a heterogeneous
            # compound.
            cluster_composition_dict = {}
            for site in bonds[sites[0]]:
                if Element(site.specie) in cluster_composition_dict:
                    cluster_composition_dict[Element(site.specie)] += 1
                else:
                    cluster_composition_dict[Element(site.specie)] = 1
            uniform = True
            if len(cluster_composition_dict):
                cmp = Composition.from_dict(cluster_composition_dict)
                if cmp.reduced_formula != structure.composition.reduced_formula:
                    uniform = False
            if not uniform:
                type = 'heterogeneous'
            else:
                # Make a 2x2x2 supercell and recalculate the
                # cluster's new size. If the new cluster size is
                # the same as the old size, it is a non-periodic
                # molecule. If it is 2x as big, it's a 1D chain.
                # If it's 4x as big, it is a layered material.
                old_cluster_size = len(bonds[sites[0]])
                structure.make_supercell(2)
                sites = structure.sites
                bonds = {}
                for site in sites:
                    bonds[site] = []

                for i in range(len(sites)):
                    site_1 = sites[i]
                    for site_2 in sites[i+1:]:
                        if (site_1.distance(site_2) <
                                float(Element(site_1.specie).atomic_radius
                                + Element(site_2.specie).atomic_radius) * 1.1):
                            bonds[site_1].append(site_2)
                            bonds[site_2].append(site_1)

                cluster_terminated = False
                while not cluster_terminated:
                    original_cluster_size = len(bonds[sites[0]])
                    for site in bonds[sites[0]]:
                        bonds[sites[0]] += [
                            s for s in bonds[site] if s not in bonds[sites[0]]]
                    if len(bonds[sites[0]]) == original_cluster_size:
                        cluster_terminated = True

                if len(bonds[sites[0]]) != 4 * old_cluster_size:
                    type = 'molecular'
                else:
                    type = 'layered'

    if write_poscar_from_cluster:
        Structure.from_sites(original_cluster).to('POSCAR', 'POSCAR')

    return type


def add_vacuum(delta, cut=0.9):
    """
    Adds vacuum to a POSCAR.

    delta = vacuum thickness in Angstroms
    cut = height above which atoms will need to be fixed. Defaults to
    0.9.
    """

    # Fix the POSCAR to put bottom atoms (even if they are above the
    # current vacuum layer) at 0.0.

    structure = Structure.from_file('POSCAR')
    n_sites = structure.num_sites
    poscar_lines = open('POSCAR').readlines()
    with open('POSCAR', 'w') as poscar:
        for line in poscar_lines[:8]:
            poscar.write(line)
        for line in poscar_lines[8:8+n_sites]:
            split_line = line.split()
            if float(split_line[2]) > cut:
                new_line = ' '.join([split_line[0], split_line[1],
                            str(float(split_line[2]) - 1.0)])
            else:
                new_line = ' '.join(split_line)
            poscar.write(new_line + '\n')

    min_z = 1
    for site in structure.sites:
        if site._fcoords[2] > cut:
            height = site._fcoords[2] - 1
        else:
            height = site._fcoords[2]
        if height < min_z:
            min_z = height

    translation = SymmOp.from_rotation_and_translation(
        translation_vec=(0, 0, -min_z))
    structure.apply_operation(translation, fractional=True)
    structure.to('POSCAR', 'POSCAR')
    with open('POSCAR', 'r') as poscar:
        poscar_lines = poscar.readlines()
    atom_lines = []
    for i in range(8, 8+n_sites):
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
#    if fixable:
#        add_factor = 1.0 - min(addables)
#    else:
    add_factor = 0.0
    new_atom_lines = []
    for atom_line in atom_lines:
        new_atom_line_2 = str(float(atom_line[2]) + add_factor)
        if float(new_atom_line_2) >= 1.0:
            new_atom_line_2 = str(float(new_atom_line_2) - 1.0)
        new_atom_lines.append('{} {} {}'.format(atom_line[0], atom_line[1],
                                                new_atom_line_2))
    with open('POSCAR', 'w') as poscar:
        for line in poscar_lines[0:8]:
            poscar.write(line)
        for new_atom_line in new_atom_lines:
            poscar.write('{}\n'.format(new_atom_line))

    # Open files and read in values from POSCAR
    old_poscar = open('POSCAR', 'r')
    new_poscar = open('new_POSCAR', 'w')
    oldlines = old_poscar.readlines()
    name = oldlines[0].split()[0]
    lattice_constant = oldlines[1].strip()
    a_lat_par = [float(x) for x in oldlines[2].split()]
    b_lat_par = [float(y) for y in oldlines[3].split()]
    c_lat_par = [abs(float(z)) for z in oldlines[4].split()]
    elements = oldlines[5].split()
    stoichiometry = oldlines[6].split()
    coordinate_type = oldlines[7].strip()

    # Elongate c-vector by delta

    save = float(c_lat_par[2])
    c_length = float(c_lat_par[2]) * float(lattice_constant)
    c_length_plus_delta = c_length + float(delta)
    c_lat_par[2] = c_length_plus_delta / float(lattice_constant)
    scalar = c_lat_par[2] / save

    # Create list of atom coordinates and adjust their z-coordinate on
    # the fly

    atoms = []
    for i in range(8, 8+n_sites):
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
    os.remove('POSCAR')

    new_lines = open('new_POSCAR').readlines()
    with open('POSCAR', 'w') as poscar:
        for line in new_lines:
            poscar.write(line)
    old_poscar.close()
    os.remove('new_POSCAR')


def write_potcar(pot_path=POTENTIAL_PATH, types='None'):
    """
    Writes a POTCAR file based on a list of types.

    types = list of same length as number of elements containing specifications
    for the kind of potential desired for each element. If no special potential
    is desired, just enter '', or leave types = 'None'.
    (['pv', '', '3'])
    """

    poscar = open('POSCAR', 'r')
    lines = poscar.readlines()
    elements = lines[5].split()
    poscar.close()

    potcar_symbols = loadfn(os.path.join(PACKAGE_PATH, 'potcar_symbols.yaml'))

    if types == 'None':
        types = [potcar_symbols[elt].replace(elt, '').replace('_', '')
                 for elt in elements]

    potentials = []
    for i in range(len(elements)):
        if types[i] == '':
            pass
        else:
            elements[i] += '_{}'.format(types[i])

        # If specified pseudopotential doesn't exist, try other variations.
        if os.path.exists('{}/{}/POTCAR'.format(pot_path, elements[i])):
            pass
        else:
            print('Potential file for {} does not exist. Looking for best'\
                  'variation... '.format(elements[i]))
            if types[i] == 'regular':
                length = 0
            else:
                length = len(types[i]) + 1
                elements[i] = elements[i][:-length]
            elements[i] += '_sv'
            if os.path.exists('{}/{}/POTCAR'.format(
                    pot_path, elements[i])):
                print('Found one! {} will work.'.format(elements[i]))
            else:
                elements[i] = elements[i][:-3]
                elements[i] += '_pv'
                if os.path.exists('{}/{}/POTCAR'.format(
                        pot_path, elements[i])):
                    print('Found one! {} will work.'.format(elements[i]))
                else:
                    elements[i] = elements[i][:-3]
                    elements[i] += '_3'
                    if os.path.exists('{}/{}/POTCAR'.format(
                            pot_path, elements[i])):
                        print('Found one! {} will work.'.format(elements[i]))
                    else:
                        elements[i] = elements[i][:-2]
                        if os.path.exists('{}/{}/POTCAR'.format(
                                pot_path, elements[i])):
                            print(('Found one! {} will '
                                   'work.'.format(elements[i])))
                        else:
                            print('No pseudopotential found'
                                   ' for {}'.format(elements[i]))

    # Create paths, open files, and write files to POTCAR for each potential.
    for element in elements:
        potentials.append('{}/{}/POTCAR'.format(pot_path, element))
    outfile = open('POTCAR', 'w')
    for potential in potentials:
        infile = open(potential)
        for line in infile:
            outfile.write(line)
        infile.close()
    outfile.close()


def write_pbs_runjob(name, nnodes, nprocessors, pmem, walltime, binary):
    """
    writes a runjob based on a name, nnodes, nprocessors, walltime, and
    binary. Designed for runjobs on the Hennig group_list on HiperGator
    1 (PBS).
    """
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
    runjob.write('mpirun {} > job.log\n\n'.format(binary))
    runjob.write('echo \'Done.\'\n')
    runjob.close()


def write_slurm_runjob(name, ntasks, pmem, walltime, binary):
    """
    writes a runjob based on a name, nnodes, nprocessors, walltime, and
    binary. Designed for runjobs on the Hennig group_list on HiperGator
    1 (PBS).
    """

    nnodes = int(np.ceil(float(ntasks) / 32.0))

    runjob = open('runjob', 'w')
    runjob.write('#!/bin/bash\n')
    runjob.write('#SBATCH --job-name={}\n'.format(name))
    runjob.write('#SBATCH -o out_%j.log\n')
    runjob.write('#SBATCH -e err_%j.log\n')
    runjob.write('#SBATCH --qos=hennig-b\n')
    runjob.write('#SBATCH --nodes={}\n'.format(nnodes))
    runjob.write('#SBATCH --ntasks={}\n'.format(ntasks))
    runjob.write('#SBATCH --mem-per-cpu={}\n'.format(pmem))
    runjob.write('#SBATCH -t {}\n\n'.format(walltime))
    runjob.write('cd $SLURM_SUBMIT_DIR\n\n')
    runjob.write('module load intel/2016.0.109\n')
    runjob.write('module load openmpi/1.10.1\n')
    runjob.write('module load vasp/5.4.1\n\n')
    runjob.write('mpirun {} > job.log\n\n'.format(binary))
    runjob.write('echo \'Done.\'\n')
    runjob.close()
