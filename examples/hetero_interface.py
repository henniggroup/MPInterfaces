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
substrate.to(fmt='poscar', filename='POSCAR_substrate_initial.vasp')
mat2d.to(fmt='poscar', filename='POSCAR_mat2d_initial.vasp')        
    
# get the matching substrate and 2D material lattices
uv_substrate, uv_mat2d = get_matching_lattices(substrate, mat2d,
                                                max_area = 200,
                                                max_mismatch = 0.01,
                                                max_angle_diff = 1,
                                                r1r2_tol = 0.02)

# map the intial slabs to the newly found matching lattices
substrate_latt = Lattice( np.array(
    [ uv_substrate[0][:],
    uv_substrate[1][:],
    substrate.lattice.matrix[2,:]
    ] ))    
mat2d_latt = Lattice(np.array(
    [ uv_mat2d[0][:],
    uv_mat2d[1][:],
    mat2d.lattice.matrix[2,:]
    ] ))
_, __, scell = substrate.lattice.find_mapping(substrate_latt,
                                              ltol = 0.01,atol = 1)
substrate.make_supercell(scell)
_, __, scell = mat2d.lattice.find_mapping(mat2d_latt,
                                            ltol = 0.01,atol = 1)
mat2d.make_supercell(scell)
substrate.to(fmt='poscar', filename='POSCAR_substrate_matching.vasp')
mat2d.to(fmt='poscar', filename='POSCAR_mat2d_matching.vasp')            

# modify the substrate lattice so that the 2d material can be
# grafted on top of it
lmap = Lattice( np.array(
    [ mat2d.lattice.matrix[0,:],
        mat2d.lattice.matrix[1,:],
    substrate.lattice.matrix[2,:]
    ] ) )    
substrate.modify_lattice(lmap)

# put the 2d material on top of the substrate
# seperation between the 2dmaterial and the substrate in angstroms    
seperation = 5
substrate_top_z = np.max(np.array([site.coords for site in substrate])[:,2])
# shift origin and vector
shift = substrate.lattice.matrix[2,:]
shift = shift/np.linalg.norm(shift) * seperation
origin = np.array([0,0, substrate_top_z])    
for site in mat2d:
        new_coords = site.coords - origin  +  shift
        substrate.append(site.specie, new_coords, coords_are_cartesian=True)
    
substrate.to(fmt='poscar', filename='POSCAR_final.vasp')        
    

#structure from materials project, use your own key
#gaas = get_struct_from_mp('GaAs', MAPI_KEY="dwvz2XCFUEI9fJiR")
#sa_gaas = SpacegroupAnalyzer(gaas)
#gaas_cvn = sa_gaas.get_conventional_standard_structure()
#gaas_cvn.to(fmt='poscar', filename='POSCAR_GaAs.vasp')

#cdte = get_struct_from_mp('CdTe', MAPI_KEY="dwvz2XCFUEI9fJiR")
#sa_cdte = SpacegroupAnalyzer(cdte)
#cdte_cvn = sa_cdte.get_conventional_standard_structure()
#cdte_cvn.to(fmt='poscar', filename='POSCAR_CdTe.vasp')    
    
