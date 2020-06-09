# coding: utf-8
# Copywrite (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Utility script to find if endpoint slabs are matching with substrate

Usage: Provide inputs <path to poscar substrate>, <path to poscar 2D> and
       match constraints using for the search

       python get_matching_matrix.py
"""
# import the lma module
from mpinterfaces import transformations
from pymatgen import Structure, Lattice
import numpy as np

#### Inputs for lattice matching #####
# substrate file path
sub_file = '../test_files/POSCAR_substrate'
# twod file path
twod_file = '../test_files/POSCAR_2D'
sub = Structure.from_file(sub_file)
twod = Structure.from_file(twod_file)

# provide input match constraints as a dictionary
match_constraints = {'max_area':90, 'max_mismatch':0.06, 'max_angle_diff':2,
                     'r1r2_tol':0.04, 'separation': 2.2, 'nlayers_substrate':1,
                     'nlayers_2d':1, 'sd_layers':0, 'best_match': 'area'}
#########################################

# Lowest area (default) or minimum uv mismatched lattice match
# is returned based on 'best_match' option in match_constraints
# best_match should be either 'area' or 'mismatch'
iface, n_sub, z_ub = transformations.run_lat_match(sub, twod, match_constraints)

if iface is None:
    print ('Current twod lattice matches with substrate. No changes needed!')
else:
    strct = Structure(iface.lattice, iface.species, iface.frac_coords)
    strct.to(filename='../test_files/POSCAR_iface', fmt='poscar')
