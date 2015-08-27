from __future__ import division, unicode_literals, print_function

"""
Compute the reduced matching lattice vectors for heterostructure
interfaces as described in the paper by Zur and McGill:
Journal of Applied Physics 55, 378 (1984); doi: 10.1063/1.333084
"""

__author__="Kiran Mathew, Arunima Singh"

import sys
import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.io.vasp.inputs import Poscar

from mpinterfaces.calibrate import CalibrateSlab
from mpinterfaces.transformations import *
from mpinterfaces.utils import *

seperation = 5 # in angstroms
nlayers_2d = 2
nlayers_substrate = 2

substrate_slab = slab_from_file([0,0,1], 'POSCAR_substrate')
mat2d_slab = slab_from_file([0,0,1], 'POSCAR_2D')
# get the in-plane lattice aligned slabs
substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    substrate_slab,
    mat2d_slab,
    max_area = 400,
    max_mismatch = 0.05,
    max_angle_diff = 1,
    r1r2_tol = 0.01)
substrate_slab_aligned.to(fmt='poscar',
                          filename='POSCAR_substrate_aligned.vasp')
mat2d_slab_aligned.to(fmt='poscar',
                      filename='POSCAR_mat2d_aligned.vasp')
# merge substrate and mat2d in all possible ways
hetero_interfaces = generate_all_configs(mat2d_slab_aligned,
                                         substrate_slab_aligned,
                                         nlayers_2d, nlayers_substrate,
                                         seperation )
# generate all poscars
for i, iface in enumerate(hetero_interfaces):
    sd_flags = CalibrateSlab.set_sd_flags2d(
        interface=iface,
        n_layers=nlayers_2d+nlayers_substrate)
    poscar = Poscar(iface, selective_dynamics=sd_flags)
    poscar.write_file(
        filename='POSCAR_final_{}.vasp'.format(i))   
    
