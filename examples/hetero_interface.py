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
                         filename='POSCAR_final_'+str(i)+'_'+str(j)+'.vasp')
            
if __name__ == '__main__':
    # initial bulk unit cells
    # mind: the hkl specification in the slab generation is wrt the
    # initial structure
    # example: GaAs used as the substrate and
    # CdTe used as the 2d material. Both conventional unit cells
    substrate_bulk = Structure.from_file('POSCAR.mp-2534_GaAs')
    mat2d_bulk = Structure.from_file('POSCAR.mp-406_CdTe')    
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
    substrate = Interface(substrate_bulk,
                          hkl = [1,1,0],
                          min_thick = 10,
                          min_vac = 25,
                          primitive = False, from_ase = True)
    mat2d = Interface(mat2d_bulk,
                      hkl = [1,1,0],
                      min_thick = 2,
                      min_vac = 0,
                      primitive = False, from_ase = True)
    substrate.to(fmt='poscar',
                 filename='POSCAR_substrate_initial.vasp')
    mat2d.to(fmt='poscar',
             filename='POSCAR_mat2d_initial.vasp')
    # get the matching substrate and 2D material lattices
    uv_substrate, uv_mat2d = get_matching_lattices(substrate, mat2d,
                                              max_area = 200,
                                              max_mismatch = 0.01,
                                              max_angle_diff = 1,
                                              r1r2_tol = 0.02)
    # map the intial slabs to the newly found matching lattices
    substrate_latt = Lattice( np.array(
                                    [
                                    uv_substrate[0][:],
                                    uv_substrate[1][:],
                                    substrate.lattice.matrix[2,:]
                                    ] ))    
    mat2d_latt = Lattice( np.array(
                                   [
                                   uv_mat2d[0][:],
                                   uv_mat2d[1][:],
                                   mat2d.lattice.matrix[2,:]
                                   ] ))
    _, __, scell = substrate.lattice.find_mapping(substrate_latt,
                                              ltol = 0.01,
                                              atol = 1)
    substrate.make_supercell(scell)
    _, __, scell = mat2d.lattice.find_mapping(mat2d_latt,
                                            ltol = 0.01,
                                            atol = 1)
    mat2d.make_supercell(scell)
    substrate.to(fmt='poscar', filename='POSCAR_substrate_matching.vasp')
    mat2d.to(fmt='poscar', filename='POSCAR_mat2d_matching.vasp')
    # modify the substrate lattice so that the 2d material can be
    # grafted on top of it
    lmap = Lattice( np.array(
        [
            mat2d.lattice.matrix[0,:],
            mat2d.lattice.matrix[1,:],
            substrate.lattice.matrix[2,:]
        ] ) )    
    substrate.modify_lattice(lmap)
    # put the 2d material on top of the substrate
    # seperation between the 2dmaterial and the substrate in angstroms
    # generate all configurations for each unique atom in the
    # specified layers of the 2d material and the substrate
    seperation = 5
    nlayers_2d = 2
    nlayers_substrate = 2
    generate_all_configs(mat2d, substrate,
                         nlayers_2d, nlayers_substrate,
                         seperation )
    #structure from materials project, use your own key
    #gaas = get_struct_from_mp('GaAs', MAPI_KEY="")
    #sa_gaas = SpacegroupAnalyzer(gaas)
    #gaas_cvn = sa_gaas.get_conventional_standard_structure()
    #gaas_cvn.to(fmt='poscar', filename='POSCAR_GaAs.vasp')
    #cdte = get_struct_from_mp('CdTe', MAPI_KEY="")
    #sa_cdte = SpacegroupAnalyzer(cdte)
    #cdte_cvn = sa_cdte.get_conventional_standard_structure()
    #cdte_cvn.to(fmt='poscar', filename='POSCAR_CdTe.vasp')    
    
