from __future__ import division, unicode_literals, print_function

"""
Compute the reduced matching lattice vectors for heterostructure
interfaces as described in the paper by Zur and McGill:
Journal of Applied Physics 55, 378 (1984); doi: 10.1063/1.333084
"""

__author__="Kiran Mathew, Arunima Singh"

import sys
from math import sqrt
from copy import copy

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import Slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from mpinterfaces import get_struct_from_mp
from mpinterfaces.interface import Interface
from mpinterfaces.transformations import *


def generate_all_configs(mat2d, substrate,
                         nlayers_2d = 2, nlayers_substrate = 2 ,
                         seperation = 5 ):
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
    if not(mat2d and substrate):
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
                                       for site in substrate])[:,2])
    mat_2d_bottom = np.min(np.array([site.coords
                                     for site in mat2d])[:,2])
    # shift normal to the surface by 'seperation'
    surface_normal = substrate.lattice.matrix[2,:]
    origin = np.array([0,0,substrate_top_z])
    shift_normal = surface_normal/np.linalg.norm(surface_normal) * seperation
    #generate all possible interfaces, one for each combination of
    # unique substrate and unique 2d materials site in the layers .i.e
    # an interface structure for each parallel shift
    #interface = 2D material + substrate
    for i, coord_i in enumerate(coords_uniq_sub):
        for j, coord_j in enumerate(coords_uniq_2d):
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
            interface.to(fmt='poscar',
                         filename='POSCAR_final_{0}_{1}.vasp'.format(i,j))

            
def get_aligned_lattices(slab_sub, slab_2d, max_area = 200,
                         max_mismatch = 0.05,
                         max_angle_diff = 1, r1r2_tol = 0.2):
    """
    given the 2 slab structures and the alignment paramters, return
    slab structures with lattices that are aligned with respect to each
    other
    """
    # get the matching substrate and 2D material lattices
    uv_substrate, uv_mat2d = get_matching_lattices(slab_sub, slab_2d,
                                              max_area = max_area,
                                              max_mismatch = max_mismatch,
                                              max_angle_diff = max_angle_diff,
                                              r1r2_tol = r1r2_tol)
    if not uv_substrate and not uv_mat2d:
        print("no matching u and v, trying adjusting the parameters")
        sys.exit()
    substrate = Structure.from_sites(slab_sub)
    mat2d = Structure.from_sites(slab_2d)
    # map the intial slabs to the newly found matching lattices
    substrate_latt = Lattice( np.array(
                                    [
                                    uv_substrate[0][:],
                                    uv_substrate[1][:],
                                    substrate.lattice.matrix[2,:]
                                    ] ))
    # to avoid numerical issues with find_mapping
    mat2d_fake_c = mat2d.lattice.matrix[2,:] / np.linalg.norm(mat2d.lattice.matrix[2,:]) * 5.0
    mat2d_latt = Lattice( np.array(
                                   [
                                   uv_mat2d[0][:],
                                   uv_mat2d[1][:],
                                   mat2d_fake_c
                                   ] ))
    mat2d_latt_fake = Lattice( np.array(
                                   [
                                   mat2d.lattice.matrix[0,:],
                                   mat2d.lattice.matrix[1,:],
                                   mat2d_fake_c
                                   ] ))    
    _, __, scell = substrate.lattice.find_mapping(substrate_latt,
                                              ltol = 0.05,
                                              atol = 1)
    scell[2] = np.array([0,0,1]) 
    substrate.make_supercell(scell)
    _, __, scell = mat2d_latt_fake.find_mapping(mat2d_latt,
                                            ltol = 0.05,
                                            atol = 1)
    scell[2] = np.array([0,0,1]) 
    mat2d.make_supercell(scell)
    # modify the substrate lattice so that the 2d material can be
    # grafted on top of it
    lmap = Lattice( np.array(
        [
            mat2d.lattice.matrix[0,:],
            mat2d.lattice.matrix[1,:],
            substrate.lattice.matrix[2,:]
        ] ) )    
    substrate.modify_lattice(lmap)
    return substrate, mat2d
    

if __name__ == '__main__':
    #
    # BULK
    #
    # mind:
    #     the hkl specification in the slab generation is wrt the
    #     initial structure
    #structure from materials project, use your own key
    substrate_bulk = get_struct_from_mp('Ag')
    sa_sub = SpacegroupAnalyzer(substrate_bulk)
    substrate_bulk = sa_sub.get_conventional_standard_structure()
    #mat2d_bulk = get_struct_from_mp('GaN')
    #sa_mat2d = SpacegroupAnalyzer(mat2d_bulk)
    #mat2d_bulk = sa_mat2d.get_conventional_standard_structure()
    #substrate_bulk.to(fmt='poscar', filename='POSCAR_substarte_bulk.vasp')    
    #mat2d_bulk.to(fmt='poscar', filename='POSCAR_mat2d_bulk.vasp')    
    #
    # SLABS
    #
    # initialize the substrate and 2d material slabs, constructed
    # from the respective bulk unit cells
    # notes:
    #      ase backend used to ensure that the generated slabs have
    #      orthogonal z axis
    #      2d material vacuum spacing = 0
    #      keep in mind that the 2d material will be put on top of
    #      the subtrate in the substrate's vacuum space.
    #      So ensure that the substrate vacuum spacing is sufficietly
    #      large enough to contain the 2d material
    # substrate of minimum thickness 10 A and 25 A minimum vacuum
    substrate_slab = Interface(substrate_bulk,
                          hkl = [1,1,1],
                          min_thick = 10,
                          min_vac = 25,
                          primitive = False, from_ase = True)
    # GaN 001 substrate structure read from file and converted to slab object
    # Note: must specify the read in slab's hkl
    mat2d_input = Structure.from_file('POSCAR_2D')
    mat2d_hkl = [0,0,1]
    mat2d_slab = Slab(mat2d_input.lattice,
                      mat2d_input.species_and_occu,
                      mat2d_input.frac_coords,
                      mat2d_hkl,
                      Structure.from_sites(mat2d_input,to_unit_cell=True),
                      shift=0,
                      scale_factor=np.eye(3, dtype=np.int),
                      site_properties=mat2d_input.site_properties)
    # create from bulk a substrate with thickness 3 A and no vacuum
    #mat2d_slab = Interface(mat2d_bulk,
    #                      hkl = [0,0,1],
    #                      min_thick = 3,
    #                      min_vac = 1,
    #                      primitive = False, from_ase = True)
    substrate_slab.to(fmt='poscar', filename='POSCAR_substrate_slab.vasp')
    mat2d_slab.to(fmt='poscar', filename='POSCAR_mat2d_slab.vasp')
    #
    # LATTICE ALIGNMENT
    #
    # get the in-plane lattice aligned slabs
    substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(substrate_slab,
                                                                      mat2d_slab,
                                                                      max_area = 100,
                                                                      max_mismatch = 0.06,
                                                                      max_angle_diff = 1,
                                                                      r1r2_tol = 0.1)
    substrate_slab_aligned.to(fmt='poscar',
                              filename='POSCAR_substrate_aligned.vasp')
    mat2d_slab_aligned.to(fmt='poscar',
                          filename='POSCAR_mat2d_aligned.vasp')
    #
    # JOINT SUPERCELL
    #
    # merge substrate and mat2d in all possible ways
    seperation = 5 # in angstroms
    nlayers_2d = 2
    nlayers_substrate = 2
    generate_all_configs(mat2d_slab_aligned, substrate_slab_aligned,
                         nlayers_2d, nlayers_substrate,
                         seperation )
