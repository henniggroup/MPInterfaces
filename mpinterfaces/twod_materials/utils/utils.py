from __future__ import print_function, division, unicode_literals

import itertools as it
import math
import os

import numpy as np

from monty.serialization import loadfn

from pymatgen.core.composition import Composition
from pymatgen.core.operations import SymmOp
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from mpinterfaces import twod_materials
from mpinterfaces.twod_materials import POTENTIAL_PATH

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__version__ = "1.6"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


def is_converged(directory):
    """
    Check if a relaxation has converged.

    Args:
        directory (str): path to directory to check.

    Returns:
        boolean. Whether or not the job is converged.
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
    Slurm systems.

    Args:
        directory (str): *absolute* path to the directory to check.

    Returns:
        string version of SLURM job status.
            'None' = No running job in this directory
    """

    os.system("squeue -o '%.18i %.9P %.16j %.8u %.8T %.10M %.9l %.6D %R %Z'"
              ">> all_jobs.txt")
    lines = open('all_jobs.txt').readlines()
    job_state = None
    for i in range(len(lines)):
        if directory in lines[i]:
            job_state = lines[i][4]
            break

    os.system('rm all_jobs.txt')

    return job_state


def get_magmom_string():
    """
    Based on a POSCAR, returns the string required for the MAGMOM
    setting in the INCAR. Initializes transition metals with 6.0
    bohr magneton and all others with 0.5.

    Returns:
        string. e.g. '1*6.0 3*0.5'
    """

    magmoms = []
    poscar_lines = open('POSCAR').readlines()
    elements = poscar_lines[5].split()
    amounts = poscar_lines[6].split()
    for i in range(len(elements)):
        if Element(elements[i]).is_transition_metal:
            magmoms.append('{}*6.0'.format(amounts[i]))
        else:
            magmoms.append('{}*0.5'.format(amounts[i]))
    return ' '.join(magmoms)


def get_spacing(filename='POSCAR', cut=0.9):
    """
    Returns the interlayer spacing for a 2D material.

    Args:
        filename (str): 'CONTCAR' or 'POSCAR', whichever file to
            check.
        cut (float): a fractional z-coordinate that must be within
            the vacuum region.

    Returns:
        float. Spacing in Angstroms.
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
    spacing = ((1.0 + min_height) - max_height) * \
              abs(float(c_axis[2])) * float(lattice_parameter[0])

    return spacing


def get_rotation_matrix(axis, theta):
    """
    Find the rotation matrix associated with counterclockwise rotation
    about the given axis by theta radians.
    Credit: http://stackoverflow.com/users/190597/unutbu

    Args:
        axis (list): rotation axis of the form [x, y, z]
        theta (float): rotational angle in radians

    Returns:
        array. Rotation matrix.
    """

    axis = np.array(list(axis))
    axis = axis / np.linalg.norm(axis)
    axis *= -np.sin(theta/2.0)
    a = np.cos(theta/2.0)
    b, c, d = tuple(axis.tolist())
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

    Args:
        structure (Structure): Pymatgen Structure object to rotate.

    Returns:
        structure. Rotated to align c-axis along [001].
    """

    c = structure.lattice._matrix[2]
    z = [0, 0, 1]
    axis = np.cross(c, z)
    if not(axis[0] == 0 and axis[1] == 0):
        theta = (np.arccos(np.dot(c, z) / (np.linalg.norm(c) * np.linalg.norm(z))))
        R = get_rotation_matrix(axis, theta)
        rotation = SymmOp.from_rotation_and_translation(rotation_matrix=R)
        structure.apply_operation(rotation)
    return structure


def get_structure_type(structure, write_poscar_from_cluster=False):
    """
    This is a topology-scaling algorithm used to describe the
    periodicity of bonded clusters in a bulk structure.

    Args:
        structure (structure): Pymatgen structure object to classify.
        write_poscar_from_cluster (bool): Set to True to write a
            POSCAR from the sites in the cluster.

    Returns:
        string. 'molecular' (0D), 'chain', 'layered', 'heterogeneous'
            (intercalated 3D), or 'conventional' (3D)
    """

    # The conventional standard structure is much easier to work
    # with.

    structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()

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

    Args:
        delta (float): vacuum thickness in Angstroms
        cut (delta): fractional height above which atoms will
            need to be fixed. Defaults to 0.9.
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
    # fixable = False
    addables = []
    for atom_line_2 in atom_line_2s:
        if float(atom_line_2) > cut or\
                (float(atom_line_2) < 0.0 and 1.0 + float(atom_line_2) > cut):
            if float(atom_line_2) < 0.0 and 1.0 + float(atom_line_2) > cut:
                atom_line_2 = float(atom_line_2) + 1.0
            addables.append(atom_line_2)
            fixable = True
    # if fixable:
    #     add_factor = 1.0 - min(addables)
    #  else:
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

    Args:
        pot_path (str): can be changed to override default location
            of POTCAR files.
        types (list): list of same length as number of elements
            containing specifications for the kind of potential
            desired for each element, e.g. ['Na_pv', 'O_s']. If
            left as 'None', uses the defaults in the
            'potcar_symbols.yaml' file in the package root.
    """

    if pot_path == '/path/to/POTCAR/files':
        # This means the config.yaml file has not been set up.
        pass
    else:
        poscar = open('POSCAR', 'r')
        lines = poscar.readlines()
        elements = lines[5].split()
        poscar.close()

        ROOT_PATH = os.path.dirname(twod_materials.__file__)
        potcar_symbols = loadfn(os.path.join(ROOT_PATH, 'potcar_symbols.yaml'))

        if types == 'None':
            sorted_types = [potcar_symbols[elt] for elt in elements]
        else:
            sorted_types = []
            for elt in elements:
                for t in types:
                    if t.split('_')[0] == elt:
                        sorted_types.append(t)

        potentials = []
        for i in range(len(elements)):
            if types[i] == '':
                pass
            else:
                elements[i] += '_{}'.format(types[i])

        # Create paths, open files, and write files to
        # POTCAR for each potential.
        for potential in sorted_types:
            potentials.append('{}/{}/POTCAR'.format(pot_path, potential))
        outfile = open('POTCAR', 'w')
        for potential in potentials:
            infile = open(potential)
            for line in infile:
                outfile.write(line)
            infile.close()
        outfile.close()


def write_circle_mesh_kpoints(center=[0, 0, 0], radius=0.1, resolution=20):
    """
    Create a circular mesh of k-points centered around a specific
    k-point and write it to the KPOINTS file. Non-circular meshes
    are not supported, but shouldn't be too hard to code. All
    k-point weights are 1.

    Args:
        center (list): x, y, and z coordinates of mesh center.
            Defaults to Gamma.
        radius (float): Size of the mesh in inverse Angstroms.
        resolution (int): Number of mesh divisions along the
            radius in the 3 primary directions.
    """

    kpoints = []
    step = radius / resolution

    for i in range(-resolution, resolution):
        for j in range(-resolution, resolution):
            if i**2 + j**2 <= resolution**2:
                kpoints.append([str(center[0]+step*i),
                                str(center[1]+step*j), '0', '1'])
    with open('KPOINTS', 'w') as kpts:
        kpts.write('KPOINTS\n{}\ndirect\n'.format(len(kpoints)))
        for kpt in kpoints:
            kpts.write(' '.join(kpt))
            kpts.write('\n')


def get_markovian_path(points):
    """
    Calculates the shortest path connecting an array of 2D
    points.

    Args:
        points (list): list/array of points of the format
            [[x_1, y_1, z_1], [x_2, y_2, z_2], ...]

    Returns:
        list: A sorted list of the points in order on the markovian path.
    """

    def dist(x, y):
        return math.hypot(y[0] - x[0], y[1] - x[1])

    paths = [p for p in it.permutations(points)]
    path_distances = [
        sum(map(lambda x: dist(x[0], x[1]), zip(p[:-1], p[1:])))
        for p in paths]
    min_index = np.argmin(path_distances)

    return paths[min_index]


def remove_z_kpoints():
    """
    Strips all linemode k-points from the KPOINTS file that include a
    z-component, since these are not relevant for 2D materials.
    """

    kpoint_lines = open('KPOINTS').readlines()

    twod_kpoints = []
    labels = {}
    i = 4

    while i < len(kpoint_lines):
        kpt_1 = kpoint_lines[i].split()
        kpt_2 = kpoint_lines[i+1].split()
        if float(kpt_1[2]) == 0.0 and [float(kpt_1[0]),
                                       float(kpt_1[1])] not in twod_kpoints:
            twod_kpoints.append([float(kpt_1[0]), float(kpt_1[1])])
            labels[kpt_1[4]] = [float(kpt_1[0]), float(kpt_1[1])]

        if float(kpt_2[2]) == 0.0 and [float(kpt_2[0]),
                                       float(kpt_2[1])] not in twod_kpoints:
            twod_kpoints.append([float(kpt_2[0]), float(kpt_2[1])])
            labels[kpt_2[4]] = [float(kpt_2[0]), float(kpt_2[1])]
        i += 3

    kpath = get_markovian_path(twod_kpoints)

    with open('KPOINTS', 'w') as kpts:
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

            kpts.write(' '.join([str(kpath[i][0]), str(kpath[i][1]), '0.0 !',
                                label_1]))
            kpts.write('\n')
            kpts.write(' '.join([str(kpt_2[0]), str(kpt_2[1]), '0.0 !',
                                label_2]))
            kpts.write('\n\n')


def write_pbs_runjob(name, nnodes, nprocessors, pmem, walltime, binary):
    """
    writes a runjob based on a name, nnodes, nprocessors, walltime,
    and binary. Designed for runjobs on the Hennig group_list on
    HiperGator 1 (PBS).

    Args:
        name (str): job name.
        nnodes (int): number of requested nodes.
        nprocessors (int): number of requested processors.
        pmem (str): requested memory including units, e.g. '1600mb'.
        walltime (str): requested wall time, hh:mm:ss e.g. '2:00:00'.
        binary (str): absolute path to binary to run.
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
    2 (SLURM).

    Args:
        name (str): job name.
        ntasks (int): total number of requested processors.
        pmem (str): requested memory including units, e.g. '1600mb'.
        walltime (str): requested wall time, hh:mm:ss e.g. '2:00:00'.
        binary (str): absolute path to binary to run.
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
