# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Compute the reduced matching lattice vectors for heterostructure
interfaces as described in the paper by Zur and McGill:
Journal of Applied Physics 55, 378 (1984); doi: 10.1063/1.333084
"""

from six.moves import range
import sys
from math import sqrt
import numpy as np

from pymatgen import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import \
    RotationTransformation

__author__ = "Kiran Mathew, Arunima Singh, V. S. Chaitanya Kolluru"
__copyright__ = "Copyright 2018, Henniggroup"
__maintainer__ = "V. S. Chaitanya Kolluru"
__email__ = "chaitanya.ismu@gmail.com"
__status__ = "Production"
__date__ = "May 19, 2020"


def get_trans_matrices(n):
    """
    Returns a list of 2x2 transformation matrices for the
    given supercell

    n: size of the supercell (area wise only)
    """

    factors = []
    for i in range(1, n + 1):
        if n % i == 0:
            factors.append(i)

    trans_matrices = []
    for i in factors:
        m = n // i
        tm = [[[i, j], [0, m]] for j in range(m)]
        trans_matrices.append(tm)
    all_tm = []
    for i in range(0, len(trans_matrices)):
        each_set = trans_matrices[i]
        for j in range(0,len(each_set)):
            all_tm.append(each_set[j])
    return all_tm


def get_uv(ab, t_mat):
    """
    Return u and v, the supercell lattice vectors obtained through the
    transformation matrix
    """
    u = np.array(ab[0]) * t_mat[0][0] + np.array(ab[1]) * t_mat[0][1]
    v = np.array(ab[1]) * t_mat[1][1]
    return [u, v]


def get_reduced_uv(uv, tm):
    """
    Returns reduced lattice vectors
    """
    is_not_reduced = True
    u = np.array(uv[0])
    v = np.array(uv[1])
    tm1 = np.array(tm)
    u1 = u.copy()
    v1 = v.copy()
    while is_not_reduced:
        if np.dot(u, v) < 0:
            v = -v
            tm1[1] = -tm1[1]
        if np.linalg.norm(u) > np.linalg.norm(v):
            u1 = v.copy()
            v1 = u.copy()
            tm1c = tm1.copy()
            tm1[0], tm1[1] = tm1c[1], tm1c[0]
        elif np.linalg.norm(v) > np.linalg.norm(u + v):
            v1 = v + u
            tm1[1] = tm1[1] + tm1[0]
        elif np.linalg.norm(v) > np.linalg.norm(u - v):
            v1 = v - u
            tm1[1] = tm1[1] - tm1[0]
        else:
            is_not_reduced = False
        u = u1.copy()
        v = v1.copy()
    return [u, v], tm1


def reduced_supercell_vectors(ab, n):
    """
    Returns all possible reduced in-plane lattice vectors and
    transition matrices for the given starting unit cell lattice
    vectors(ab) and the supercell size n
    """
    uv_list = []
    tm_list = []
    for r_tm in get_trans_matrices(n):
        uv = get_uv(ab, r_tm)
        uv_r, tm0 = get_reduced_uv(uv, r_tm)
        uv_list.append(uv_r)
        tm_list.append(tm0)
    return uv_list, tm_list


def get_r_list(area1, area2, max_area, tol=0.02):
    """
    Returns a list of r1 and r2 values that satisfies:
    r1/r2 = area2/area1 with the constraints:
    r1 <= Area_max/area1 and r2 <= Area_max/area2
    r1 and r2 corresponds to the supercell sizes of the 2 interfaces
    that align them
    """
    r_list = []
    rmax1 = int(max_area / area1)
    rmax2 = int(max_area / area2)
    print('rmax1, rmax2: {0}, {1}\n'.format(rmax1, rmax2))
    for r1 in range(1, rmax1 + 1):
        for r2 in range(1, rmax2 + 1):
            if abs(float(r1) * area1 - float(r2) * area2) / max_area <= tol:
                r_list.append([r1, r2])
    return r_list


def get_mismatch(a, b):
    """
    Percentage mistmatch between the lattice vectors a and b
    """
    a = np.array(a)
    b = np.array(b)
    return np.linalg.norm(b) / np.linalg.norm(a) - 1


def get_angle(a, b):
    """
    Angle between lattice vectors a and b in degrees
    """
    a = np.array(a)
    b = np.array(b)
    return np.arccos(
        np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b)) * 180 / np.pi


def surface_area(cell):
    """
    Calculates the surface area of the Cell
    """
    m = cell.lattice.matrix
    return np.linalg.norm(np.cross(m[0], m[1]))


def get_area(uv):
    """
    Returns area of the parallelogram, given a and b
    """
    a = uv[0]
    b = uv[1]
    return np.linalg.norm(np.cross(a, b))


def remove_duplicates(uv_list, tm_list):
    """
    Remove duplicates based on a, b, alpha matching.
    """
    new_uv_list = []
    new_tm_list = []
    for sup_lat_n1, tm_n1 in zip(uv_list, tm_list):
        a1          = [np.linalg.norm(i[0]) for i in sup_lat_n1]
        b1          = [np.linalg.norm(i[1]) for i in sup_lat_n1]
        angles      = [get_angle(i[0], i[1]) for i in sup_lat_n1]
        n1_lattices = [(a, b, alpha) for a, b, alpha in zip(a1, b1, angles)]
        for lat in n1_lattices:
            zround      = np.array(n1_lattices).round(1)
            zlist       = zround.tolist()
            zstr        = np.array([str(j) for j in zlist])
            zu, zind    = np.unique(zstr, return_index = True)
            unq_sup_lat = [sup_lat_n1[i] for i in zind]
            unq_tm      = [tm_n1[i] for i in zind]
        new_uv_list.append(unq_sup_lat)
        new_tm_list.append(unq_tm)
    return new_uv_list, new_tm_list


def get_matching_lattices(iface1, iface2, max_area=100,
                          max_mismatch=0.01, max_angle_diff=1,
                          r1r2_tol=0.02, opt=False, best_match='area'):
## This function gives: "return uv_opt[0], uv_opt[1]"

    """
    computes a list of matching reduced lattice vectors that satify
    the max_area, max_mismatch and max_anglele_diff criteria
    """
    if iface1 is None and iface2 is None:
        # test : the numbers from the paper
        a1 = 5.653
        a2 = 6.481
        # for 100 plane
        ab1 = [[0, a1 / 2, -a1 / 2], [0, a1 / 2, a1 / 2]]
        ab2 = [[0, a2 / 2, -a2 / 2], [0, a2 / 2, a2 / 2]]
        area1 = a1 ** 2 / 2
        area2 = a2 ** 2 / 2
        # for 110 plane
        ab1 = [[a1 / 2, -a1 / 2, 0], [0, 0, a1]]
        ab2 = [[a2 / 2, -a2 / 2, 0], [0, 0, a2]]
        area1 = a1 ** 2 / sqrt(2)
        area2 = a2 ** 2 / sqrt(2)
        # for 111 surface
        # ab1 = [ [a1/2, 0, a1/2], [a1/2, a1/2, 0]]
        # ab2 = [ [a2/2, 0, a2/2], [a2/2, a2/2, 0]]
        # area1 = a1**2 * sqrt(3)/4 #/ 2 /sqrt(2)
        # area2 = a2**2 * sqrt(3)/4 #/ 2 / sqrt(2)
    else:
        area1 = surface_area(iface1)
        area2 = surface_area(iface2)
        # ab1 is list of two lattice vectors that define the substrate lattice
        # ab2 is that of 2d material
        ab1 = [iface1.lattice.matrix[0, :], iface1.lattice.matrix[1, :]]
        ab2 = [iface2.lattice.matrix[0, :], iface2.lattice.matrix[1, :]]
    #print('initial values:\nuv1:\n{0}\nuv2:\n{1}\n '.format(ab1, ab2))
    r_list = get_r_list(area1, area2, max_area, tol=r1r2_tol)
    if not r_list:
        print('r_list is empty. Try increasing the max surface '
                'area or/and the other tolerance parameters')
        return None, None
        #sys.exit()
    found = []
    #print('searching ...')
    uv1_list, tm1_list, uv2_list, tm2_list = [], [], [], []
    for r1r2 in r_list:
        x1, y1 = reduced_supercell_vectors(ab1, r1r2[0])
        uv1_list.append(x1)
        tm1_list.append(y1)
        x2, y2 = reduced_supercell_vectors(ab2, r1r2[1])
        uv2_list.append(x2)
        tm2_list.append(y2)
        if not uv1_list and not uv2_list:
            continue
    new_uv1, new_tm1 = remove_duplicates(uv1_list, tm1_list)
    new_uv2, new_tm2 = remove_duplicates(uv2_list, tm2_list)
    for sup_lat_n1, sup_lat_n2 in zip(new_uv1, new_uv2):
        for i, uv1 in enumerate(sup_lat_n1):
            for j, uv2 in enumerate(sup_lat_n2):
                u_mismatch = get_mismatch(uv1[0], uv2[0])
                v_mismatch = get_mismatch(uv1[1], uv2[1])
                angle1 = get_angle(uv1[0], uv1[1])
                angle2 = get_angle(uv2[0], uv2[1])
                angle_mismatch = abs(angle1 - angle2)
                area1 = get_area(uv1)
                area2 = get_area(uv2)
                if abs(u_mismatch) < max_mismatch and abs(
                        v_mismatch) < max_mismatch:
                    if angle_mismatch < max_angle_diff:
                        found.append((uv1, uv2, max(area1, area2), u_mismatch,
                                      v_mismatch, angle_mismatch))
        # to increase the speed of the algorithm
        # Since the algorithm searches by increasing order from r_list,
        # lowest area match is found first
        # Or if past 59 seconds from the search loop , exit
                if found:
                    break   # stop searching when first match is found

    if found:
        print('\nMATCH FOUND\n')
        if best_match == 'area': # sort based on area and return lowest area
            uv_all = sorted(found, key=lambda x: x[2])
        elif best_match == 'mismatch': # sort based on average of uv mismatches
            uv_all = sorted(found, key=lambda x: (abs(x[3]) + abs(x[4])) / 2)
        uv_opt = uv_all[0]                # min. area match
        print('Best match:\nuv1:\n{0}\nuv2:\n{1}\narea:\n{2}\n'.format(
            uv_opt[0], uv_opt[1], uv_opt[2]))
        print('Lattice mismatch[u, v & alpha]:\n{0} \%, {1} \%, \
        {2} degrees\n'.format(uv_opt[3]*100, uv_opt[4]*100, uv_opt[5]))

        #print('\nSmallest area matched uv\n')
        return uv_opt[0], uv_opt[1]
    else:
        print('\n NO MATCH FOUND\n')
        return None, None


def get_uniq_layercoords(struct, nlayers, top=True):
    """
    returns the coordinates of unique sites in the top or bottom
    nlayers of the given structure.
    Args:
        struct: input structure
        nlayers: number of layers
        top: top or bottom layers, default is top layer
    Return:
        numpy array of unique coordinates
    """
    coords = np.array([site.coords for site in struct])
    z = coords[:, 2]
    z = np.around(z, decimals=4)
    zu, zuind = np.unique(z, return_index=True)
    if top:
        z_nthlayer = z[zuind[-nlayers]]
        zfilter = (z >= z_nthlayer)
    else:
        z_nthlayer = z[zuind[nlayers - 1]]
        zfilter = (z <= z_nthlayer)
    # site indices in the layers
    indices_layers = np.argwhere(zfilter).ravel()
    sa = SpacegroupAnalyzer(struct)
    symm_data = sa.get_symmetry_dataset()
    # equivalency mapping for the structure
    # i'th site in the struct equivalent to eq_struct[i]'th site
    eq_struct = symm_data["equivalent_atoms"]
    # equivalency mapping for the layers
    eq_layers = eq_struct[indices_layers]
    # site indices of unique atoms in the layers
    __, ueq_layers_indices = np.unique(eq_layers, return_index=True)
    # print(ueq_layers_indices)
    indices_uniq = indices_layers[ueq_layers_indices]
    # coordinates of the unique atoms in the layers
    return coords[indices_uniq]

def get_interface(substrate, mat2d, nlayers_2d=2, nlayers_substrate=2,
                                 separation=5):
    """
    For the given lattice matched 2D material and substrate structures,
    this functions computes all unique sites in the interface layers
    and subsequently generates all possible unique 2d/substrate
    interfaces and writes the corresponding poscar files
    Args:
        mat2d: Lattice and symmetry-matched 2D material structure
        substrate: Lattice and symmetry-matched 2D substrate structure
        nlayers_substrate: number of substrate layers
        nlayers_2d: number of 2d material layers
        separation: separation between the substrate and the 2d
                    material
    Returns:
        None
    TODO: give additional random placement of 2D material on substrate
    """
    # immediate exit if no structures
    if not (mat2d and substrate):
        print("no structures. aborting lattice match ...")
        return None
        #sys.exit()
    # unique site coordinates in the substrate top layers
    coords_uniq_sub = get_uniq_layercoords(substrate,
                                           nlayers_substrate,
                                           top=True)
    # unique site coordinates in the 2D material bottom layers
    coords_uniq_2d = get_uniq_layercoords(mat2d,
                                          nlayers_2d,
                                          top=False)
    substrate_top_z = np.max(np.array([site.coords
                                       for site in substrate])[:, 2])
    mat_2d_bottom = np.min(np.array([site.coords
                                     for site in mat2d])[:, 2])
    # shift normal to the surface by 'seperation'
    surface_normal = substrate.lattice.matrix[2, :]
    origin = np.array([0, 0, substrate_top_z])
    shift_normal = surface_normal / np.linalg.norm(surface_normal) * separation
    # generate all possible interfaces, one for each combination of
    # unique substrate and unique 2d materials site in the layers .i.e
    # an interface structure for each parallel shift
    # interface = 2D material + substrate
    interface = substrate.copy()
    shift_parallel = coords_uniq_sub[0] - coords_uniq_2d[0]
    shift_parallel[2] = 0
    shift_net = shift_normal - shift_parallel

    # generate new coords for 2D material to be added to substrate
    new_coords = []
    inds = []
    mat_species = []
    for ind, site in enumerate(mat2d):
        new_coord = site.coords
        new_coord[2] = site.coords[2] - mat_2d_bottom
        new_coord = new_coord + origin + shift_net
        new_coords.append(new_coord)
        inds.append(ind)
        mat_species.append(site.specie)

    inds = np.array(inds) + len(substrate)
    # insert mat2d coords and species at the end of the interface using index
    # Not using Structure.append method as it seems to disrupt atoms order
    for i, specie, coord in zip(inds, mat_species, new_coords):
        interface.insert(i, specie, coord, coords_are_cartesian=True)

    return interface


def get_aligned_lattices(slab_sub, slab_2d, max_area=200,
                         max_mismatch=0.05,
                         max_angle_diff=1, r1r2_tol=0.2, best_match='area'):
    """
    given the 2 slab structures and the alignment paramters, return
    slab structures with lattices that are aligned with respect to each
    other
    """

    # get the matching substrate and 2D material lattices
    uv_substrate, uv_mat2d = get_matching_lattices(
                                                slab_sub, slab_2d,
                                                max_area=max_area,
                                                max_mismatch=max_mismatch,
                                                max_angle_diff=max_angle_diff,
                                                r1r2_tol=r1r2_tol,
                                                best_match=best_match)
    if not uv_substrate and not uv_mat2d:
        print("no matching u and v, trying adjusting the parameters")
        return None, None
        #sys.exit()

    substrate = Structure.from_sites(slab_sub)
    mat2d = Structure.from_sites(slab_2d)

    # map the intial slabs to the newly found matching lattices
    substrate_latt = Lattice(np.array(
            [
                uv_substrate[0][:],
                uv_substrate[1][:],
                substrate.lattice.matrix[2, :]
            ]))
    # to avoid numerical issues with find_mapping
    mat2d_fake_c = mat2d.lattice.matrix[2, :] / np.linalg.norm(
            mat2d.lattice.matrix[2, :]) * 5.0
    mat2d_latt = Lattice(np.array(
            [
                uv_mat2d[0][:],
                uv_mat2d[1][:],
                mat2d_fake_c
            ]))

    mat2d_latt_fake = Lattice(np.array(
            [
                mat2d.lattice.matrix[0, :],
                mat2d.lattice.matrix[1, :],
                mat2d_fake_c
            ]))

    # Supercell matrix for primitive lattices -> match lattices
    _, __, scell_1 = substrate.lattice.find_mapping(substrate_latt,
                                                      ltol=0.05,
                                                      atol=max_angle_diff)
    _, __, scell_2 = mat2d_latt_fake.find_mapping(mat2d_latt,
                                                    ltol=0.05,
                                                    atol=max_angle_diff)
    scell_1[2] = np.array([0, 0, 1])
    substrate.make_supercell(scell_1)
    scell_2[2] = np.array([0, 0, 1])
    mat2d.make_supercell(scell_2)

    # modify the substrate lattice
    lmap = Lattice(np.array(
            [
                substrate.lattice.matrix[0, :],
                substrate.lattice.matrix[1, :],
                mat2d.lattice.matrix[2, :]
            ]))
    mat2d.lattice = lmap

    return substrate, mat2d

def rotate_to_principal_directions(cell):
    """
    Author: Benjamin Revard

    Rotates the cell into the principal directions. That is, lattice vector
    a is parallel to the Cartesian x-axis, lattice vector b lies in the
    Cartesian x-y plane and the z-component of lattice vector c is
    positive.
    Note: this method doesn't change the fractional coordinates of the
    sites. However, the Cartesian coordinates may be changed.
    """

    # rotate about the z-axis to align a vertically with the x-axis
    rotation = RotationTransformation(
        [0, 0, 1], 180 - (180/np.pi)*np.arctan2(
                cell.lattice.matrix[0][1],
                cell.lattice.matrix[0][0]))
    new_structure = rotation.apply_transformation(cell)
    cell.lattice = new_structure.lattice

    # rotate about the y-axis to make a parallel to the x-axis
    rotation = RotationTransformation(
            [0, 1, 0], (180/np.pi)*np.arctan2(
                cell.lattice.matrix[0][2],
                cell.lattice.matrix[0][0]))
    new_structure = rotation.apply_transformation(cell)
    cell.lattice = new_structure.lattice
    # rotate about the x-axis to make b lie in the x-y plane
    rotation = RotationTransformation(
            [1, 0, 0], 180 - (180/np.pi)*np.arctan2(
                cell.lattice.matrix[1][2],
                cell.lattice.matrix[1][1]))
    new_structure = rotation.apply_transformation(cell)
    cell.lattice = new_structure.lattice

    # make sure they are all pointing in positive directions
    if cell.lattice.matrix[0][0] < 0:
        # rotate about y-axis to make a positive
        rotation = RotationTransformation([0, 1, 0], 180)
        new_structure = rotation.apply_transformation(cell)
        cell.lattice = new_structure.lattice
    if cell.lattice.matrix[1][1] < 0:
        # rotate about x-axis to make b positive
        rotation = RotationTransformation([1, 0, 0], 180)
        new_structure = rotation.apply_transformation(cell)
        cell.lattice = new_structure.lattice
    if cell.lattice.matrix[2][2] < 0:
        # mirror c across the x-y plane to make it positive
        # a and b
        a = cell.lattice.matrix[0]
        b = cell.lattice.matrix[1]
        # the components of c
        cx = cell.lattice.matrix[2][0]
        cy = cell.lattice.matrix[2][1]
        cz = -1*cell.lattice.matrix[2][2]
        cell.lattice = Lattice([a, b, [cx, cy, cz]])

def run_lat_match(substrate, twod_layer, match_constraints):

    '''
    Runs the lattice matching algorithm on a substrate and a 2D materials

    Args:

        substrate - substrate Cell

        twod_layer - 2D layer Cell

        match_constraints - dictionary containing max area, max mismatch of u
        or v, max angle difference, area ratio tolerence, seperation at the
        interface, number of layers for substrate and 2D. Ex: {'max_area':200,
        'max_mismatch':0.05, 'max_angle_diff':2, 'r1r2_tol':0.06, 'separation':
        3, 'nlayers_substrate':1, 'nlayers_2d':1, 'sd_layers':1}
    '''
    # variables from the keys
    max_area = match_constraints['max_area']
    max_mismatch = match_constraints['max_mismatch']
    max_angle_diff = match_constraints['max_angle_diff']
    r1r2_tol = match_constraints['r1r2_tol']
    separation = match_constraints['separation']
    nlayers_substrate = match_constraints['nlayers_substrate']
    nlayers_2d = match_constraints['nlayers_2d']
    sd_layers = match_constraints['sd_layers']
    best_match = match_constraints['best_match']

    twod_prim = twod_layer.get_primitive_structure()
    substrate_prim = substrate.get_primitive_structure()
    n_prim_sub = substrate_prim.num_sites

    try:
        #get aligned lattices
        sub, mat2d = get_aligned_lattices(
                            substrate_prim,
                            twod_prim,
                            max_area=max_area,
                            max_mismatch=max_mismatch,
                            max_angle_diff=max_angle_diff,
                            r1r2_tol=r1r2_tol,
                            best_match=best_match)
        rotate_to_principal_directions(sub)
        rotate_to_principal_directions(mat2d)
        # sorts atoms wrt electronegativity
        # use this order in POTCAR
        sub.sort()
        mat2d.sort()
        n_aligned_sub = sub.num_sites
        scell_size = n_aligned_sub / n_prim_sub
    except:
        print ('Lattice match failed due to singular matrix generation for supercell')
        return None, None, None

    #merge substrate and mat2d in all possible ways
    hetero_interfaces = None
    if sub and mat2d:
        try:
            hetero_interface = get_interface(sub, mat2d,
                                     nlayers_2d, nlayers_substrate,
                                     separation)
        except:
            print('Lattice match failed at get_interface')
            return None, None, None

        z_coords_sub = sub.frac_coords[:, 2]
        z_unique, z_inds = np.unique(z_coords_sub, return_index=True)
        if sd_layers == 0: # freeze all substrate atoms
            sd_index = n_aligned_sub - 1
        else:    # relax top layer of substrate atoms
            sd_index = z_inds[len(z_inds)-sd_layers] - 1

    if hetero_interface:
        return  hetero_interface, n_aligned_sub, sd_index
    else:
        return None, None, None
