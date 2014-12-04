"""
TO DO: WRITE IN FUNCTION FORMAT TAKING XYZ FILE AS INPUT --- OBJECTIVE NOW TO CONSTRUCT A POSCAR OF A DI-MOLECULE , DIMER, ETC. FEATURES TO ADD SITES AS WELL, SOME CLEANUP REQUIRED WITH INTERFACE

LONG TERM: testing to create a new molecule from existing molecules, can be taken forward to any hybrid of two structures

 METHOD: first trying to make lead di-acetate from acetic acid: take single xyz of acetic acid, convert to molecule1 and molecule2,

  apply translation and rotation to molecule2 , convert back to structure? to remove hydrogen sites, add lead and join/adsorb together through nearest neighbor atom 
	# DESIRED USER INPUTS: Bridge atom index, rotation and translation to original molecule on the copy
"""
import sys
from pymatgen.core.structure import Structure, Molecule, Lattice
import numpy as np
from pymatgen.core.operations import SymmOp

## IGNORE CAN BE WORKED ON LATER , skip to do : work with reflection , converting to rotation and translate for now - ATTEMPT TO USE REFLECTION 
##ref_site_mol1 = mol1.cart_coords[i]
##origin_mirror = ref_site_mol1 + mirror_shift
##normal_vec = origin_mirror - ref_site_mol1
##reflect_shift = SymmOp.reflection(normal_vec, origin_mirror)
##translate_vec = np.array([0, -7, 0])#np.array([shift[0], shift[1], shift[2]])
##print translate_vec

def combine_mols(mol1, mol2, shift, angle_rot):
    array_shift =  np.random.rand(1,3) + 1.0 #np.array(shift)
    #re-center on bridge and rotate about bridge atom
    #first translation to recenter at bridge index
    Ox = angle_rot[0]
    Oy = angle_rot[1]
    Oz = angle_rot[2]

    distances = [mol1.get_distance(0,i) for i in range(len(mol1.cart_coords))]
    farthest = np.argmax(distances)
    #print 'farthest from 0', farthest
    #vector pointing from 0th index atom to the farthest atom of mol1
    vec = mol1.cart_coords[farthest] - mol1.cart_coords[0]
    #rotation
    #rotation_x = trans_and_rot.from_axis_angle_and_translation((1, 0, 0), Ox)
    #rotation_y = trans_and_rot.from_axis_angle_and_translation((0, 1, 0), Oy)
    #rotation_z = trans_and_rot.from_axis_angle_and_translation((0, 0, 1), Oz)
    #mol2.apply_operation(rotation_x)
    #mol2.apply_operation(rotation_y)
    #mol2.apply_operation(rotation_z)
    #print "\n After rotation coordinates: \n", mol.cart_coords
    op = SymmOp.from_origin_axis_angle(mol2.cart_coords[0], axis=[-vec[0], vec[1], vec[2]], angle=90)
    mol2.apply_operation(op)
    new_coords = mol1.cart_coords + array_shift #- mol2.cart_coords[0]
    #for i, coord in enumerate(mol1.cart_coords):
    #    if i == 0:
    #        coord = coord + array_shift
    #    new_coords.append(coord)
    mol1 =  mol1.__class__(mol1.species_and_occu, new_coords,
                          charge=mol1._charge,
                          spin_multiplicity=mol1._spin_multiplicity,
                          site_properties=mol1.site_properties)
    ## CAN BE IGNORED, ATTEMPT TO USE REFLECTION - continued 
    ##mol2.apply_operation(reflect_shift)
    #final shift translation
    #length = len(mol2.cart_coords)
    #index = [] #for index input needed for translate sites
    #for x in xrange(length):
    #    index.append(x) 
    #mol2.translate_sites(index, array_shift)
    #print mol2
    #remove sites from molecule , removing hydrogens located at site (8-1 = 7)
    mol2.remove_sites([7])
    mol1.remove_sites([7])
    combine_mol_sites = mol1.sites + mol2.sites
    #combine molecule units
    combine_mol = Molecule.from_sites(combine_mol_sites, validate_proximity=True) 

    return combine_mol


if __name__=='__main__':
    shift = [0.0, -5.0, 0.0]
    angle_rot =[180, 0, 0]    
    mol1 = Molecule.from_file("acetic_acid.xyz")
    mol2 = Molecule.from_file("acetic_acid.xyz")
    #print mol1
    #print mol2

    combine_mol = combine_mols(mol1, mol2, shift, angle_rot)
    combine_mol.to('xyz', 'combo.xyz')
    print combine_mol

    # get boxed structure of molecule in a predefined default box 
    combine_mol_struct = combine_mol.get_boxed_structure(13, 13, 13)  
    print combine_mol_struct

    #TEST SITE FOR Pb atom , not final 
    combine_mol_struct.append('Pb', [0.417020, 0.2, 0.5]) 
    combine_mol_struct.to(fmt= "poscar", filename= "POSCAR_diacetate_boxed.vasp")

  





                
                
