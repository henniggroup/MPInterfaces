# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
create nanoparticle using wulff construction
"""

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from mpinterfaces import get_struct_from_mp
from mpinterfaces.nanoparticle import Nanoparticle

# -----------------------------------
# nanoparticle specifications
# -----------------------------------
# max radius in angstroms
rmax = 15
# surface families to be chopped off
hkl_family = [(1, 0, 0), (1, 1, 1)]
# surfac energies could be in any units, will be normalized
surface_energies = [28, 25]

# -----------------------------------
# initial structure
# -----------------------------------
# caution: set the structure wrt which the the miller indices are
# specified. use your own key
structure = get_struct_from_mp('PbS', MAPI_KEY="")
# primitive ---> conventional cell
sa = SpacegroupAnalyzer(structure)
structure_conventional = sa.get_conventional_standard_structure()

# -----------------------------------
# create nanoparticle
# -----------------------------------
nanoparticle = Nanoparticle(structure_conventional, rmax=rmax,
                            hkl_family=hkl_family,
                            surface_energies=surface_energies)
nanoparticle.create()
nanoparticle.to(fmt='xyz', filename='nanoparticle.xyz')
