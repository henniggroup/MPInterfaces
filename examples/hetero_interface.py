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

__author__ = "Kiran Mathew, Arunima Singh"

from mpinterfaces.calibrate import CalibrateSlab
from mpinterfaces.interface import Interface
from mpinterfaces.transformations import *
from mpinterfaces.utils import *

separation = 3  # in angstroms
nlayers_2d = 2
nlayers_substrate = 2

substrate_bulk = Structure.from_file('POSCAR_substrate')
# substrate_bulk = get_struct_from_mp('Ag')
sa_sub = SpacegroupAnalyzer(substrate_bulk)
substrate_bulk = sa_sub.get_conventional_standard_structure()
substrate_slab = Interface(substrate_bulk,
                           hkl=[1, 1, 1],
                           min_thick=10,
                           min_vac=25,
                           primitive=False, from_ase=True)
# substrate_slab = slab_from_file([0,0,1], 'POSCAR_substrate')
mat2d_slab = slab_from_file([0, 0, 1], 'POSCAR_2D')
# get the in-plane lattice aligned slabs
# substrate_slab.to(fmt='poscar', filename='POSCAR_substrate_slab.vasp')
mat2d_slab.to(fmt='poscar', filename='POSCAR_mat2d_slab.vasp')
# selective dynamics flag
sd_flags = CalibrateSlab.set_sd_flags(
    interface=substrate_slab,
    n_layers=nlayers_substrate,
    top=True, bottom=False)
poscar = Poscar(substrate_slab, selective_dynamics=sd_flags)
poscar.write_file(filename='POSCAR_substrate_slab.vasp')
# get aligned lattices
substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
    substrate_slab,
    mat2d_slab,
    max_area=400,
    max_mismatch=0.05,
    max_angle_diff=1,
    r1r2_tol=0.01)
substrate_slab_aligned.to(fmt='poscar',
                          filename='POSCAR_substrate_aligned.vasp')
mat2d_slab_aligned.to(fmt='poscar',
                      filename='POSCAR_mat2d_aligned.vasp')
# merge substrate and mat2d in all possible ways
hetero_interfaces = generate_all_configs(mat2d_slab_aligned,
                                         substrate_slab_aligned,
                                         nlayers_2d, nlayers_substrate,
                                         separation)
# generate all poscars
for i, iface in enumerate(hetero_interfaces):
    sd_flags = CalibrateSlab.set_sd_flags(
        interface=iface,
        n_layers=nlayers_2d + nlayers_substrate,
        top=True, bottom=False)
    poscar = Poscar(iface, selective_dynamics=sd_flags)
    poscar.write_file(
        filename='POSCAR_final_{}.vasp'.format(i))
