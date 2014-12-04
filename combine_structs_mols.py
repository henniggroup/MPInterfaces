## ""TO DO: WRITE IN FUNCTION FORMAT TAKING XYZ FILE AS INPUT --- OBJECTIVE NOW TO CONSTRUCT A POSCAR OF A DI-MOLECULE , DIMER, ETC. FEATURES TO ADD SITES AS WELL, SOME CLEANUP REQUIRED WITH INTERFACE 
##LONG TERM: testing to create a new molecule from existing molecules, can be taken forward to any hybrid of two structures 
 # METHOD: first trying to make lead di-acetate from acetic acid: take single xyz of acetic acid, convert to molecule1 and molecule2,
 # apply translation and rotation to molecule2 , convert back to structure? to remove hydrogen sites, add lead and join/adsorb together through nearest neighbor atom 
	# DESIRED USER INPUTS: Bridge atom index, rotation and translation to original molecule on the copy  ""


# In[63]:

from pymatgen.core.structure import Structure, Molecule, Lattice
import numpy as np
from pymatgen.core.operations import SymmOp as trans_and_rot


# In[64]:

shift = [0.0, -5.0, 0.0]
angle_rot =[180, 0, 0]
i = 0
array_shift = np.array(shift)
#print mirror_shift
print array_shift


# In[65]:

mol1 = Molecule.from_file("/home/josh/Research/POSCARs/xyzs/acetic_acid_new.xyz")
mol2 = Molecule.from_file("/home/josh/Research/POSCARs/xyzs/acetic_acid_new.xyz")
print mol1


# In[66]:

## IGNORE CAN BE WORKED ON LATER , skip to do : work with reflection , converting to rotation and translate for now - ATTEMPT TO USE REFLECTION 
##ref_site_mol1 = mol1.cart_coords[i]
##origin_mirror = ref_site_mol1 + mirror_shift
##normal_vec = origin_mirror - ref_site_mol1


# In[67]:

##reflect_shift = SymmOp.reflection(normal_vec, origin_mirror)
##translate_vec = np.array([0, -7, 0])#np.array([shift[0], shift[1], shift[2]])
##print translate_vec


# In[68]:

#re-center on bridge and rotate about bridge atom
#first translation to recenter at bridge index
Ox = angle_rot[0]
Oy = angle_rot[1]
Oz = angle_rot[2]
origin = [0, 0, 0] #default origin
mol_orig_coords = mol2.cart_coords #for processing
translate_vec = origin - mol_orig_coords[i]
length = len(mol2.cart_coords)
index = [] #for index input needed for translate sites
for x in xrange(length):
    index.append(x) 
mol2.translate_sites(index, translate_vec)
#print "\n After translation coordinates: \n", mol.cart_coords
#rotation
rotation_x = trans_and_rot.from_axis_angle_and_translation((1, 0, 0), Ox)
rotation_y = trans_and_rot.from_axis_angle_and_translation((0, 1, 0), Oy)
rotation_z = trans_and_rot.from_axis_angle_and_translation((0, 0, 1), Oz)
mol2.apply_operation(rotation_x)
mol2.apply_operation(rotation_y)
mol2.apply_operation(rotation_z)
#print "\n After rotation coordinates: \n", mol.cart_coords


# In[69]:

## CAN BE IGNORED, ATTEMPT TO USE REFLECTION - continued 

##mol2.apply_operation(reflect_shift)

#final shift translation
#length = len(mol2.cart_coords)
#index = [] #for index input needed for translate sites
#for x in xrange(length):
#    index.append(x) 
#mol2.translate_sites(index, array_shift)
#print mol2


# In[70]:

#remove sites from molecule , removing hydrogens located at site (8-1 = 7)
mol2.remove_sites([7])
mol1.remove_sites([7])


# In[71]:

combine_mol_sites = mol1.sites + mol2.sites


# In[72]:

combine_mol = Molecule.from_sites(combine_mol_sites, validate_proximity=True) #combine molecule units
print combine_mol
combine_mol_struct = combine_mol.get_boxed_structure(13, 13, 13)  # get boxed structure of molecule in a predefined default box 
print combine_mol_struct
combine_mol_struct.append('Pb', [0.417020, 0.2, 0.5]) #TEST SITE FOR Pb atom , not final 
combine_mol_struct.to(fmt= "poscar", filename= "/home/josh/Research/POSCARs/POSCAR_diacetate_boxed.vasp")

  





                
                
