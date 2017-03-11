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

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Kiran Mathew, Arunima Singh"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


def get_trans_matrices(n):
    """
    yields a list of 2x2 transformation matrices for the
    given supercell
    n: size
    """

    def factors(n0):
        for i in range(1, n0 + 1):
            if n0 % i == 0:
                yield i

    for i in factors(n):
        m = n // i
        yield [[[i, j], [0, m]] for j in range(m)]


def get_uv(ab, t_mat):
    """
    return u and v, the supercell lattice vectors obtained through the
    transformation matrix
    """
    u = np.array(ab[0]) * t_mat[0][0] + np.array(ab[1]) * t_mat[0][1]
    v = np.array(ab[1]) * t_mat[1][1]
    return [u, v]


def get_reduced_uv(uv, tm):
    """
    returns reduced lattice vectors
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
    returns all possible reduced in-plane lattice vectors and
    transition matrices for the given starting unit cell lattice
    vectors(ab) and the supercell size n
    """
    uv_list = []
    tm_list = []
    for r_tm in get_trans_matrices(n):
        for tm in r_tm:
            uv = get_uv(ab, tm)
            uv, tm1 = get_reduced_uv(uv, tm)
            uv_list.append(uv)
            tm_list.append(tm1)
    return uv_list, tm_list


def get_r_list(area1, area2, max_area, tol=0.02):
    """
    returns a list of r1 and r2 values that satisfies:
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
    percentage mistmatch between the lattice vectors a and b
    """
    a = np.array(a)
    b = np.array(b)
    return np.linalg.norm(b) / np.linalg.norm(a) - 1


def get_angle(a, b):
    """
    angle between lattice vectors a and b in degrees
    """
    a = np.array(a)
    b = np.array(b)
    return np.arccos(
        np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b)) * 180 / np.pi


def get_area(uv):
    """
    area of the parallelogram
    """
    a = uv[0]
    b = uv[1]
    return np.linalg.norm(np.cross(a, b))


def get_matching_lattices(iface1, iface2, max_area=100,
                          max_mismatch=0.01, max_angle_diff=1,
                          r1r2_tol=0.02):
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
        area1 = iface1.surface_area
        area2 = iface2.surface_area
        # a, b vectors that define the surface
        ab1 = [iface1.lattice.matrix[0, :], iface1.lattice.matrix[1, :]]
        ab2 = [iface2.lattice.matrix[0, :], iface2.lattice.matrix[1, :]]
    print('initial values:\nuv1:\n{0}\nuv2:\n{1}\n '.format(ab1, ab2))
    r_list = get_r_list(area1, area2, max_area, tol=r1r2_tol)
    if not r_list:
        print(
            'r_list is empty. Try increasing the max surface area or/and the other tolerance paramaters')
        sys.exit()
    found = []
    print('searching ...')
    for r1r2 in r_list:
        uv1_list, tm1_list = reduced_supercell_vectors(ab1, r1r2[0])
        uv2_list, tm2_list = reduced_supercell_vectors(ab2, r1r2[1])
        if not uv1_list and not uv2_list:
            continue
        for i, uv1 in enumerate(uv1_list):
            for j, uv2 in enumerate(uv2_list):
                u_mismatch = get_mismatch(uv1[0], uv2[0])
                v_mismatch = get_mismatch(uv1[1], uv2[1])
                angle1 = get_angle(uv1[0], uv1[1])
                angle2 = get_angle(uv2[0], uv2[1])
                angle_mismatch = abs(angle1 - angle2)
                area1 = get_area(uv1)
                area2 = get_area(uv2)
                if abs(u_mismatch) < max_mismatch and abs(
                        v_mismatch) < max_mismatch:
                    max_angle = max(angle1, angle2)
                    min_angle = min(angle1, angle2)
                    mod_angle = max_angle % min_angle
                    is_angle_factor = False
                    if abs(mod_angle) < 0.001 or abs(
                            mod_angle - min_angle) < 0.001:
                        is_angle_factor = True
                    if angle_mismatch < max_angle_diff or is_angle_factor:
                        if angle_mismatch > max_angle_diff:
                            if angle1 > angle2:
                                uv1[1] = uv1[0] + uv1[1]
                                tm1_list[i][1] = tm1_list[i][0] + tm1_list[i][
                                    1]
                            else:
                                uv2[1] = uv2[0] + uv2[1]
                                tm2_list[j][1] = tm2_list[j][0] + tm2_list[j][
                                    1]
                        found.append((uv1, uv2, min(area1, area2), u_mismatch,
                                      v_mismatch, angle_mismatch, tm1_list[i],
                                      tm2_list[j]))
    if found:
        print('\nMATCH FOUND\n')
        uv_opt = sorted(found, key=lambda x: x[2])[0]
        print('optimum values:\nuv1:\n{0}\nuv2:\n{1}\narea:\n{2}\n'.format(
            uv_opt[0], uv_opt[1], uv_opt[2]))
        print('optimum transition matrices:\ntm1:\n{0}\ntm2:\n{1}\n'.format(
            uv_opt[6], uv_opt[7]))
        print('u,v & angle mismatches:\n{0}, {1}, {2}\n'.format(uv_opt[3],
                                                                uv_opt[4],
                                                                uv_opt[5]))
        return uv_opt[0], uv_opt[1]
    else:
        print('\n NO MATCH FOUND')
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
    z_nthlayer = z[zuind[-nlayers]]
    zfilter = (z >= z_nthlayer)
    if not top:
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


def generate_all_configs(mat2d, substrate,
                         nlayers_2d=2, nlayers_substrate=2,
                         seperation=5):
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
        seperation: seperation between the substrate and the 2d
                    material
    Returns:
        None

    TODO: give additional random placement of 2D material on substrate
    """
    # immediate exit if no structures
    if not (mat2d and substrate):
        print("no structures. aborting ...")
        sys.exit()
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
    shift_normal = surface_normal / np.linalg.norm(surface_normal) * seperation
    # generate all possible interfaces, one for each combination of
    # unique substrate and unique 2d materials site in the layers .i.e
    # an interface structure for each parallel shift
    # interface = 2D material + substrate
    hetero_interfaces = []
    for coord_i in coords_uniq_sub:
        for coord_j in coords_uniq_2d:
            interface = substrate.copy()
            shift_parallel = coord_i - coord_j
            shift_parallel[2] = 0
            shift_net = shift_normal - shift_parallel
            for site in mat2d:
                new_coords = site.coords
                new_coords[2] = site.coords[2] - mat_2d_bottom
                new_coords = new_coords + origin + shift_net
                interface.append(site.specie, new_coords,
                                 coords_are_cartesian=True)
            hetero_interfaces.append(interface)
    return hetero_interfaces


def get_aligned_lattices(slab_sub, slab_2d, max_area=200,
                         max_mismatch=0.05,
                         max_angle_diff=1, r1r2_tol=0.2):
    """
    given the 2 slab structures and the alignment paramters, return
    slab structures with lattices that are aligned with respect to each
    other
    """
    # get the matching substrate and 2D material lattices
    uv_substrate, uv_mat2d = get_matching_lattices(slab_sub, slab_2d,
                                                   max_area=max_area,
                                                   max_mismatch=max_mismatch,
                                                   max_angle_diff=max_angle_diff,
                                                   r1r2_tol=r1r2_tol)
    if not uv_substrate and not uv_mat2d:
        print("no matching u and v, trying adjusting the parameters")
        sys.exit()
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
    _, __, scell = substrate.lattice.find_mapping(substrate_latt,
                                                  ltol=0.05,
                                                  atol=1)
    scell[2] = np.array([0, 0, 1])
    substrate.make_supercell(scell)
    _, __, scell = mat2d_latt_fake.find_mapping(mat2d_latt,
                                                ltol=0.05,
                                                atol=1)
    scell[2] = np.array([0, 0, 1])
    mat2d.make_supercell(scell)
    # modify the substrate lattice so that the 2d material can be
    # grafted on top of it
    lmap = Lattice(np.array(
        [
            substrate.lattice.matrix[0, :],
            substrate.lattice.matrix[1, :],
            mat2d.lattice.matrix[2, :]
        ]))
    mat2d.modify_lattice(lmap)
    return substrate, mat2d
