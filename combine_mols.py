"""
TO DO: WRITE IN FUNCTION FORMAT TAKING XYZ FILE AS INPUT --- OBJECTIVE NOW TO CONSTRUCT A POSCAR OF A DI-MOLECULE , DIMER, ETC. FEATURES TO ADD SITES AS WELL, SOME CLEANUP REQUIRED WITH INTERFACE

LONG TERM: testing to create a new molecule from existing molecules, can be taken forward to any hybrid of two structures

 METHOD: first trying to make lead di-acetate from acetic acid: take single xyz of acetic acid, convert to molecule1 and molecule2,

  apply translation and rotation to molecule2 , convert back to structure? to remove hydrogen sites, add lead and join/adsorb together through nearest neighbor atom 
	# DESIRED USER INPUTS: Bridge atom index, rotation and translation to original molecule on the copy
"""
import sys
import copy
from pymatgen.core.structure import Structure, Molecule, Lattice
import numpy as np
from pymatgen.core.operations import SymmOp

## IGNORE CAN BE WORKED ON LATER , skip to do : work with reflection , converting to rotation and translate for now - ATTEMPT TO USE REFLECTION 
##ref_site_mol0 = mol0.cart_coords[i]
##origin_mirror = ref_site_mol0 + mirror_shift
##normal_vec = origin_mirror - ref_site_mol1
##reflect_shift = SymmOp.reflection(normal_vec, origin_mirror)
##translate_vec = np.array([0, -7, 0])#np.array([shift[0], shift[1], shift[2]])
##print translate_vec
#link={'0':[], '1':[], '2':[[6,2]]}
def combine_mols(mol0, mol1, mol2, angle, link={}, remove=[]):
    array_shift =  np.random.rand(1,3) + 1.0 
    distances = [mol0.get_distance(0,i) for i in range(len(mol0.cart_coords))]
    farthest = np.argmax(distances)
    #print 'farthest from 0', farthest
    #vector pointing from 0th index atom to the farthest atom of mol0
    vec = mol0.cart_coords[farthest] - mol0.cart_coords[0]
    op = SymmOp.from_origin_axis_angle(mol1.cart_coords[0], axis=[-vec[0], vec[1], vec[2]], angle=angle)
    mol1.apply_operation(op)
    #shift the coordiantes of mol0
    new_coords = mol0.cart_coords + array_shift #- mol1.cart_coords[0]
    mol0 =  mol0.__class__(mol0.species_and_occu, new_coords,
                          charge=mol0._charge,
                          spin_multiplicity=mol0._spin_multiplicity,
                          site_properties=mol0.site_properties)
    
    new_coords = copy.deepcopy(mol2.cart_coords)
    #connect mol2 to mol0 and mol1
    if link:
        for k, v in link.items():
            if v:
                for ind, conn in enumerate(v):
                    print 'ind, conn', ind, conn
                    coord = (mol0.cart_coords[conn[0]] + mol1.cart_coords[conn[1]] )/2
                    new_coords[ind, 0] = coord[0]
                    new_coords[ind, 1] = coord[1]
                    new_coords[ind, 2] = coord[2]                                        
                    

    mol2 =  mol2.__class__(mol2.species_and_occu, new_coords,
                          charge=mol2._charge,
                          spin_multiplicity=mol2._spin_multiplicity,
                          site_properties=mol2.site_properties)
        
    mol1.remove_sites(remove)
    mol0.remove_sites(remove)
    combine_mol_sites = mol0.sites + mol1.sites + mol2.sites
    #combine molecule units
    combine_mol = Molecule.from_sites(combine_mol_sites, validate_proximity=True) 

    return combine_mol


if __name__=='__main__':
    shift = [0.0, -5.0, 0.0]
    angle = 90 #[180, 0, 0]    
    mol0 = Molecule.from_file("acetic_acid.xyz")
    mol1 = Molecule.from_file("acetic_acid.xyz")
    mol2 = Molecule(["Pb"], [[0,0,0]])
    remove = [7]
    #print mol2
    #print mol0
    #print mol1

    combine_mol = combine_mols(mol0, mol1, mol2, angle, link={'0':[], '1':[], '2':[[6,2]]}, remove=remove)
    combine_mol.to('xyz', 'combo.xyz')
    print combine_mol

    # get boxed structure of molecule in a predefined default box 
    combine_mol_struct = combine_mol.get_boxed_structure(13, 13, 13)  
    print combine_mol_struct

    #TEST SITE FOR Pb atom , not final 
    #combine_mol_struct.append('Pb', [0.417020, 0.2, 0.5]) 
    combine_mol_struct.to(fmt= "poscar", filename= "POSCAR_diacetate_boxed.vasp")

  





                
                
