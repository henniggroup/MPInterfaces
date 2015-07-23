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

def get_all_config_pos( mat2d,substrate ):

	"""
	Inputs: Lattice-matched and symmetry-matched 2D material and 
	substrate POSCAR files
	2D Materials lattice constants adjusted to the substrate lattice constants
	Output: Returns all possible non-equivalent configurations of 2D material 
	on substrate by identifying the non-equivalent positions near the interface
	using spglib
	Gives additional random placement of 2D material on substrate (to implement)
	"""

	if mat2d is None and substrate is None:
		pos_2D=Structure.from_file('POSCAR_2D')
		pos_sub=Structure.from_file('POSCAR_substrate')
	else:
		pos_2d=mat2d
		pos_sub=substrate

	"""
	Identify the non equivalent atomic positions of the 2 POSCARs
	Extract the non equivalent atomic positions of top 2 layers of substrate
	and bottom 2 layers of 2D material (curently only top 1 layer)
	"""
	finder2d = SpacegroupAnalyzer(pos_2d)
	finderbulk= SpacegroupAnalyzer(pos_sub) 
	non_equi_2d=finder2d.get_symmetry_dataset()
	non_equi_sub=finderbulk.get_symmetry_dataset()
	set_2d=non_equi_2d["equivalent_atoms"]
	set_subs=non_equi_sub["equivalent_atoms"]

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
	#shift the atomic positions of the 2D material as all combinations of
	#substrate and 2D materials resultant atomic locations

	#coordinates of the substrate atoms
	
	coor_substrate=(np.array([site.coords for site in pos_sub]))
	
	#z coordinate array of the substrate poscar
	
	z=coor_substrate[:,2]
	
	#array of indices of the topmost layer of the substrate poscar
	
	indices_top=np.where(z==max(z))
	indices_top=indices_top[0]

	#equivalent positions list for the top most layers of atoms in substrate
	
	short_set_subs=[set_subs[i] for i in indices_top]

	#all unique atoms in the top most layer of the substrate
	
	u,indices_unique=np.unique(short_set_subs,return_index=True)
	indices_true=[indices_top[i] for i in indices_unique]

	# coordinates of the unique atoms in the top most layer
	
	coor_unique_sub=coor_substrate[indices_true]

	#coordinates of 2D material atoms
	
	coor_2d=(np.array([site.coords for site in pos_2d]))

	#z coordinate array of the 2D material poscar
	
	z=coor_2d[:,2]

	#array of indices of the bottom most layer of the substrate poscar
	
	indices_bottom=np.where(z==min(z))
	indices_bottom=indices_bottom[0]

	#equivalent positions list for the top most layers of atoms in substrate
	
	short_set_2d=[set_2d[i] for i in indices_bottom]

	#all unique atoms in the top most layer of the substrate
	
	u,indices_unique=np.unique(short_set_2d,return_index=True)
	indices_true=[indices_bottom[i] for i in indices_unique]

	# coordinates of the unique atoms in the top most layer
	
	coor_unique_2d=coor_2d[indices_true]

	#create the POSCARs of 2D material on the substrate
	
	for i in range(0,coor_unique_sub.ndim-1) :
		for j in range(0,coor_unique_2d.ndim-1) :
			shift_position=coor_unique_sub[i]-coor_unique_2d[j]
			shift_position[2]=0
			substrate_top_z = np.max(np.array([site.coords for site in substrate])[:,2])
			mat_2d_bottom=np.min(np.array([site.coords for site in mat2d])[:,2])
			#shift origin and vector
			shift = substrate.lattice.matrix[2,:]
			shift = shift/np.linalg.norm(shift) * seperation
			origin = np.array([0,0,substrate_top_z])
			print origin
			print shift
			for site in mat2d:
				new_coords=site.coords
				new_coords[2]=site.coords[2]-mat_2d_bottom
				new_coords = new_coords + origin  +  shift #-shift_position
				substrate.append(site.specie, new_coords, coords_are_cartesian=True)
        		substrate.to(fmt='poscar', filename='POSCAR_final_'+str(i)+'_'+str(j)+'.vasp')
			return None


"""
An example
"""
pos_2D=Structure.from_file('POSCAR_2D')
pos_sub=Structure.from_file('POSCAR_substrate')
a=get_all_config_pos(pos_2D,pos_sub)
