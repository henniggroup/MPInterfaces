#from __future__ import divison, unicode_literals, print_function

"""
Compute all configurations of 2D materials on substrates by identifying the Wyckoff positions of the atoms
To be used after determining the super-lattices from the hetero_inteface.py script
"""

import sys 
from math import sqrt
from copy import copy

import numpy as np

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

from mpinterfaces import get_struct_from_mp
from mpinterfaces.interface import Interface
from mpinterfaces.transformations import reduced_supercell_vectors

def get_all_config_pos(mat2d,substrate)

"""
Inputs: Lattice-matched and symmetry-matched 2D material and 
substrate POSCAR files
2D Materials lattice constants adjusted to the substrate lattice constants
Output: Returns all possible non-equivalent configurations of 2D material 
on substrate by identifying the non-equivalent positions near the interface
using spglib
Gives additional random placement of 2D material on substrate
"""

if mat2d is None and substrate is None:
	pos_2D=Structure.from_file('POSCAR_2D')
	pos_sub=Structure.from_file('POSCAR_substrate')
else:
	pos_2D=mat2d
	pos_sub=substrate

"""
Identify the non equivalent atomic positions of the 2 POSCARs
Extract the non equivalent atomic positions of top 2 layers of substrate
and bottom 2 layers of 2D material
"""

finder2D = SpacegroupAnalyzer(pos_2D)
finderbulk= SpacegroupAnalyzer(pos_sub) 
non_equi_2D=finder2D.get_symmetry_dataset()
non_equi_sub=finderbulk.get_symmetry_dataset()
2d_set=non_equi_2D["equivalent_atoms"]
subs_set=non_equi_sub["equivalent_atoms"]

"""
Align the POSCAR of 2D material on substrate with all combinations of 
non-equilavent atomic locations directly upon each other along with 
a random placement of 2D material on substrate
specify seperation of the 2D material and substrate in Angstrom
"""
seperation=5

#assuming that the POSCARs are created with z positions in ascending order
#identify array number of top two layers of substrate
# and bottom 2 layers of 2D material
# for each materials, check "equivalent_atoms" array number and delete 
#dupplicates
#save the atomic positions for each material in 2 seperate arrays
#in a loop shift the atomic positions of the 2D material as all combinations of
# substrate and 2D materials resultant atomic locations

#array index of top two layers of the substrate

B = numpy.zeros(3, int)
coor_substrate=(numpy.array([site.coords for site in pos_sub]))
for i in xrange(3):
	idx = numpy.argmax(coor_substrate)
	
