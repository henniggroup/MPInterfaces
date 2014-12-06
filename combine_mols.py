"""

read in molecules from xyz files
combine them and put it in abox and write the resulting structure to poscar file

"""

import sys
import copy
from pymatgen.core.structure import Structure, Molecule, Lattice
import numpy as np
from pymatgen.core.operations import SymmOp

#normal to 2 vectors
def get_perp_vec(vec1, vec2):
    #if the vectors are parllel, then perp_vec = (0, -z, y)
    if np.abs(np.dot(vec1, vec2) - np.linalg.norm(vec1)**2 ) < 1e-6:
        perp_vec = np.array([0, -vec1[2], vec1[1]])
    else:
        perp_vec = np.cross(vec1, vec2)
    return perp_vec                        

    
#returns distnace matrix for the given molecule
def get_distance_matrix(mol):
    nsites = len(mol.sites)
    return np.array([mol.get_distance(i,j) for i in range(nsites) for j in range(nsites)]).reshape(nsites, nsites)


#combine mols
def combine_mols(mols, cm_dist, angle={}, link={}, remove=[]):
    #get the start and end indices to define the vector that defines the molecule
    vec_indices = []
    for mol in mols:
        nsites = len(mol.sites)
        d_mat = get_distance_matrix(mol)
        max_dist = np.max(d_mat.reshape(nsites*nsites, 1))
        temp = []
        for i in range(nsites):
            if i not in temp:
                [temp.append([i,j]) for j in range(nsites) if np.abs(max_dist-d_mat[i,j]) < 1e-6]
        vec_indices.append(temp[0])

    print vec_indices

    #vectors pointing from start index atom to the farthest atom
    mol_vecs = []
    for mol,vind in enumerate(vec_indices):
        mol_vecs.append(mols[mol].cart_coords[vind[1]] - mols[mol].cart_coords[vind[0]])
    
    #move center of masses
    new_mol = mols[0]
    mov_vec = np.array([1,0,0])
    for i in range(len(mols)-1):
        cm1 = new_mol.center_of_mass
        new_cm = new_mol.center_of_mass        
        cm2 = mols[i+1].center_of_mass
        new_cm = new_cm + cm_dist[i] * mov_vec #+ np.random.rand(1,3)
        mov_vec = get_perp_vec(mol_vecs[i], mov_vec)
        mov_vec = mov_vec / np.linalg.norm(mov_vec)
        new_coords = mols[i+1].cart_coords + new_cm
        mols[i+1] = mols[i+1].__class__(mols[i+1].species_and_occu, new_coords,
                          charge=mols[i+1]._charge,
                          spin_multiplicity=mols[i+1]._spin_multiplicity,
                          site_properties=mols[i+1].site_properties)

        new_mol = Molecule.from_sites(mols[i].sites + mols[i+1].sites, validate_proximity=True) 

    #rotate the molecules around an axis that is perpendicular to the molecular axes
    if angle:
        for mol in range(len(mols)):
            for ind_key, rot in angle[str(mol)].items():
                    #print 'mol, ind_key, rot ', mol, ind_key, rot
                    perp_vec = np.cross(mol_vecs[int(ind_key)], mol_vecs[mol])
                    #if the vectors are parllel, then perp_vec = (-y, x, 0)
                    if np.abs(np.dot(mol_vecs[int(ind_key)], mol_vecs[mol]) - np.linalg.norm(mol_vecs[mol])**2 ) < 1e-6:
                        perp_vec = np.array([-mol_vecs[mol][1], mol_vecs[mol][0], 0])
                    org_pt = vec_indices[mol][0]
                    op = SymmOp.from_origin_axis_angle(mols[mol].cart_coords[org_pt], axis=perp_vec, angle=rot)
                    mols[mol].apply_operation(op)

    
    #connect the molecules together
    #in this example: connect mol2 to mol0 and mol1, the coordinates of mol2 is changed 
    new_coords = np.array([0,0,0]) #copy.deepcopy(mols[0].cart_coords)
    shift = np.array([0,0,0])
    
    if link:
        for mol in range(len(mols)):
            new_coords = copy.deepcopy(mols[mol].cart_coords)
            if link[str(mol)]:
                for ind_key, conn in link[str(mol)].items():
                    ind = int(ind_key)
                    print 'connection list for atom of index ',ind ,' of molecule ', mol, ' : ', conn
                    coord = np.array([0,0,0])
                    #if connecting the molecule mol to only one atom of jus tone another molecule
                    #then move the atom close to the atom in mol and shift the whole molecule too
                    non_neg = np.extract(np.array(conn)>0, conn)
                    if len(non_neg) == 1 and len(link[str(mol)]) == 1:
                        for j,k in enumerate(conn):
                            coord = mols[j].cart_coords[non_neg[0]] + np.random.rand(1,3) + 1.0
                        shift = coord - mols[mol].cart_coords[ind]
                    #connect the specified atoms of mol to the atoms of other molecules in the list
                    #connection means putting the atomm  of the mol at a position that is the average of the position of the atoms of the molecules given in the list
                    else:
                        for j,k in enumerate(conn):
                            if k>=0:
                                coord = coord + mols[j].cart_coords[k] #+ mols[1].cart_coords[conn[1]] )/2
                        coord = coord / len(conn)
                    new_coords[ind, 0] = coord[0]
                    new_coords[ind, 1] = coord[1]
                    new_coords[ind, 2] = coord[2]
            
                new_coords = new_coords + shift
                mols[mol] =  mols[mol].__class__(mols[mol].species_and_occu, new_coords,
                          charge=mols[mol]._charge,
                          spin_multiplicity=mols[mol]._spin_multiplicity,
                          site_properties=mols[mol].site_properties)
    
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


#test
if __name__=='__main__':
    mol0 = Molecule.from_file("acetic_acid.xyz")
    mol1 = Molecule.from_file("acetic_acid.xyz")
    mol2 = Molecule(["Pb"], [[0,0,0]])
    mols = [mol0, mol1, mol2]

    #center of mass distances in angstrom
    #example: 3 molecules and cm_dist = [4,2],
    #center of mass of mol1 is moved from mol0 in 1,0,0 direction by 4 A
    #mol2 is moved from the center of mass of the combined mol0+mol1 molecule by 2 A
    # in a direction that is perpendicular to the first moving direction and the
    #molecule vector of one of the molecules
    # for n molecules the size of cm_dist must be n-1
    cm_dist = [4, 2]

    #optional parmater
    #example: angle={'0':{}, '1':{'0':90}, '2':{} }
    #rotate mol1 with respect to mol0 by 90 degreeen around and axis that is normal
    # to the plane containing the molecule vectors of mol0 and mol1
    angle={'0':{}, '1':{'0':90}, '2':{} }
    
    
    #optional paramter
    #a dictionary describing the connection between the molecules, used if the
    #relative movemnet of the molecules throught the center of mass shift is not enough
    #the key is the index for the molecules
    #the value for each key is a list of lists wiht each list indicating how the atom of index
    # corresponding to the list index in the molecule coresponding to key  should be connectd to the atoms
    #in the list
    #if not connecting to any atom of the molecule set that index for that molecule to -1
    #example:- link = {'0':{}, '1':{}, '2':{'0':[6,2, -1]} }
    #the 0th atom of the third molecule,'2', is connected to the 6th and 2nd atoms of molecules
    #'0' and '1' respectively and no coonection to the molecule'2' (indicated by -1)
    #if not connecting a single atom of one molecule to another atom of one another molecule, connection basically means
    # putting the atom at a halfway distance between other atoms
    link = {'0':{}, '1':{}, '2':{'0':[6,2, -1]} }
        
    #list of indices of atoms to be removed from each molecule
    remove = [[7],[7],[]]
    
    #combined_mol = combine_mols(mols, cm_dist, angle, link=link, remove=remove)
    combined_mol = combine_mols(mols, cm_dist, angle=angle, link={}, remove=remove)
    combined_mol.to('xyz', 'combo.xyz')
    print combined_mol

    # get boxed structure of molecule in a predefined default box 
    combine_mol_struct = combined_mol.get_boxed_structure(13, 13, 13)  
    print combine_mol_struct

    combine_mol_struct.to(fmt= "poscar", filename= "POSCAR_diacetate_boxed.vasp")

  





                
                
