"""
TO DO: WRITE IN FUNCTION FORMAT TAKING XYZ FILE AS INPUT --- OBJECTIVE NOW TO CONSTRUCT A POSCAR OF A DI-MOLECULE , DIMER, ETC. FEATURES TO ADD SITES AS WELL, SOME CLEANUP REQUIRED WITH INTERFACE

LONG TERM: testing to create a new molecule from existing molecules, can be taken forward to any hybrid of two structures

"""

import sys
import copy
from pymatgen.core.structure import Structure, Molecule, Lattice
import numpy as np
from pymatgen.core.operations import SymmOp

def combine_mols(mols, angle, link={}, remove=[]):
    distances = [mols[0].get_distance(0,i) for i in range(len(mols[0].cart_coords))]
    farthest = np.argmax(distances)

    #vector pointing from 0th index atom to the farthest atom of mol0
    vec = mols[0].cart_coords[farthest] - mols[0].cart_coords[0]
    op = SymmOp.from_origin_axis_angle(mols[1].cart_coords[0], axis=[-vec[0], vec[1], vec[2]], angle=angle)
    mols[1].apply_operation(op)
    
    #shift the coordiantes of mol0 by some random amount
    new_coords = mols[0].cart_coords + np.random.rand(1,3) + 1.0 
    mols[0] =  mols[0].__class__(mols[0].species_and_occu, new_coords,
                          charge=mols[0]._charge,
                          spin_multiplicity=mols[0]._spin_multiplicity,
                          site_properties=mols[0].site_properties)

    
    #connect the molecules together
    #in this example: connect mol2 to mol0 and mol1
    new_coords = copy.deepcopy(mols[0].cart_coords)    
    if link:
        for i in range(len(mols)):
        #for k, v in link.items():
            if link[str(i)]:
                new_coords = copy.deepcopy(mols[2].cart_coords)
                for ind, conn in enumerate(link[str(i)]):
                    print 'ind, conn', ind, conn
                    coord = np.array([0,0,0])
                    for j,k in enumerate(conn):
                        coord = coord + mols[j].cart_coords[k] #+ mols[1].cart_coords[conn[1]] )/2
                    coord = coord / len(conn)
                    new_coords[ind, 0] = coord[0]
                    new_coords[ind, 1] = coord[1]
                    new_coords[ind, 2] = coord[2]                                        
                mols[i] =  mols[i].__class__(mols[i].species_and_occu, new_coords,
                          charge=mols[i]._charge,
                          spin_multiplicity=mols[i]._spin_multiplicity,
                          site_properties=mols[i].site_properties)
    
    #remove atoms from the molecules
    for i in range(len(mols)):
        if remove[i]:
            mols[i].remove_sites(remove[i])

    #combine the sites of all the molecules
    combine_mol_sites = mols[0].sites
    for j in range(1, len(mols)):
        combine_mol_sites = combine_mol_sites + mols[j].sites
        
    #create molecule from the combined sites
    combine_mol = Molecule.from_sites(combine_mol_sites, validate_proximity=True) 

    return combine_mol


if __name__=='__main__':
    mol0 = Molecule.from_file("acetic_acid.xyz")
    mol1 = Molecule.from_file("acetic_acid.xyz")
    mol2 = Molecule(["Pb"], [[0,0,0]])
    mols = [mol0, mol1, mol2]
    #a dictionary describing the connection between the molecules
    #the key is the index for the molecules
    #the value for each key is a list of lists wiht each list indicating how the atom of index
    # corresponding to the list index in the molecule coresponding to key  should be connectd to the atoms
    #in the list
    #example:- the 0th atom of the third molecule,'2', is connected to the 6th and 2nd atoms of molecules
    #'0' and '1' respectively
    #connection is basically putting the atom at a halfway distance between other atoms
    link = {'0':[], '1':[], '2':[[6,2]]}
    #list of indices of atoms to be removed from each molecule
    remove = [[7],[7],[]]
    #angle of rotation of mol1 with respect to mol0
    angle = 90 

    combined_mol = combine_mols(mols, angle, link=link, remove=remove)
    combined_mol.to('xyz', 'combo.xyz')
    print combined_mol

    # get boxed structure of molecule in a predefined default box 
    combine_mol_struct = combined_mol.get_boxed_structure(13, 13, 13)  
    print combine_mol_struct

    combine_mol_struct.to(fmt= "poscar", filename= "POSCAR_diacetate_boxed.vasp")

  





                
                
