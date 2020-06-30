from __future__ import print_function, division, unicode_literals

import copy
import numpy as np
import math
import itertools
import operator


from collections import Counter
from sympy import Point3D, Plane,Line3D


from pymatgen import Structure, Element, Composition
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import _pt_data

__author__ = "Joshua Paul"
__copyright__ = "Copyright 2020, Henniggroup"
__maintainer__ = "Joshua Paul"
__email__ = "joshua.thomas.paul@gmail.com"
__status__ = "Production"
__date__ = "June 6, 2020"


ELEMENT_RADII = {i: Element(i).atomic_radius for i in _pt_data}
new = {}
for i in ELEMENT_RADII:
    if ELEMENT_RADII[i]!=None:
        new[i]=ELEMENT_RADII[i]
    else:
        new[i]=-1
ELEMENT_RADII = new

def getDistMat(structure, tol):
    """
    Used to create a binary matrix of bonded (1) and 
    unbonded (0) between all atoms in the system

    inputs
    --------
        structure (Structure):  pymatgen structure object for which
                                to generate the binary matrix

        tol (float):    The tolerance for bonding. .1 == 10% increase
                        in atomic radii when looking for bonds
            
    returns
    --------
        binary_matrix (list, NxN):   NxN array of bonded sites (1) and 
                                and unbonded sites (0). 
                                N == number of sites in "structure"
    """
    
    if tol>1:
        print('WARNING: Tolerance input is greater than 100%')
        
    # Find interatomic distances)
    s = structure
    distance_matrix = np.array(s.distance_matrix)
    
    
    # Ensure diagonal will be identified as not bonded
    np.fill_diagonal(distance_matrix,100)

    radii = [ELEMENT_RADII[site.species_string] for site in s.sites]
    radiiS = []
    for r in radii:
        radiiS.append(r*(1+tol))
    radii_transpose = np.array(radii)[np.newaxis].T
    
    
    # Create matrix for sum of radii of each species, scaled by 'tol'
    radii_matrix = radiiS + radii_transpose*(1+tol)
    
    
    # Take difference between atomic distances and summed radii
    temp = distance_matrix-radii_matrix
    
    
    # Turn bonded indexes to '1' and non-bonded to '0'
    binary_matrix = (temp < 0).astype(int)
    return(binary_matrix)

def buildNetwork(binary_matrix,seed_index):
    """
    Used top build the list of atomic sites which are bonded in
    a crystal.

    inputs
    --------
        binary_matrix (list, NxN):  A matrix of 0's and 1's which
                                    represents atoms as unbonded 
                                    and bonded, respectively. Can
                                    be obtained from the "getDistMat" 
                                    function
        seed_index (int): The atomic site which is used as the starting 
                          point for building the atomic network. Typically
                          does not impact the results of the TSA, but will
                          if the system is a bipartide or has mixed 
                          dimensionality

    returns
    --------
        cluster (list, 1xN):    A list of the atomic sites which are
                                bonded in the provided "binary_matrix"
                                varible. N == number of sites in the 
                                network

    """
    
    
    # Get starting point for network
    seed = set(np.where(binary_matrix[seed_index]==1)[0])
    cluster = set(seed)
    NEW = set(seed)
    cluster.update(np.array([seed_index]))
    NEW.update(np.array([seed_index]))
    while True:
        temp_set = set()
        for n in NEW:
            # temp_set will have all atoms, without duplicates,
            # that are connected to all atoms in NEW.
            temp_set.update(np.where(binary_matrix[n]==1)[0])
        if temp_set.issubset(cluster):
            # if temp_set has no new atoms, the search is done.
            break
        else:
            NEW = temp_set - cluster # List of newly discovered atoms
            cluster.update(temp_set) # cluster is updated with new atoms
    return(cluster)

def getDim(scale,supercell):
    
    """
    Helper function for to determine periodicity of a crystal, used in 
    various functions

    inputs
    --------
        scale (float):      The ratio of final and initial network size
        supercell (int):    The supercell size used to generate the new
                            atomic network
        

    returns
    --------
        motiif (str):   The structural motiif of the crystal

    """
    
    # Check for standard scaling
    motiif_dict = {1:'molecular',supercell:'chains',\
                 supercell**2:'layered', supercell**3:'conventional'}
    if scale in motiif_dict:
        return(motiif_dict[scale])
    
    # If the structure is some intermediate, determine
    # which intermediate
    
    else:
        if scale < 1:
            motiif = 'shrunk molecular'
        elif scale < supercell:
            motiif = "mol-chain"
        elif scale < supercell**2:
            motiif = "chain-2D"
        elif scale < supercell**3:
            motiif = "2D-conv"
        else:
            motiif = 'Network size increased'
        return(motiif)

def getStructureType(entry, seed_index=0, supercell = 2,returnSS=False):
    
    """
    This is a topology-scaling algorithm used to describe the
    periodicity of bonded clusters in a bulk structure. It is 
    described in Ashton et al. in 2017. If used, please cite:
        https://doi.org/10.1103/PhysRevLett.118.106101
    
    inputs
    --------
        entry (list): A set of components necessary for the TSA.
                      Makes it easier to parallelize with this as
                      the input
                      --structure (Structure): pymatgen Structure object
                      --tol (float): The scaling for the atomic bonds
                      --mp_id (str): The label for the entry, commonly
                                     the MaterialsProject ID
        seed_index (int):   The site to use as the starting point for the
                            TSA. Typically does not impact the results, but
                            will if the structure is a bipartide or has
                            mixed dimensionality
                            
        supercell (int):    The size supercell to use for building the 
                            bonded atomic network. Typically does not 
                            matter, except in the case of a bipartide
                            crystal. In such cases, even or odd sized
                            supercells will result in the TSA identifying 
                            a clear dimensionality or "mixed" dimensionality,
                            not respectivel.
                            
                            
    returns
    --------
        type0 (str): The dimensionality of the structure
        mp_id (str): Same mp_id as entered in "entry"
        tol (float): Same float as entered in "entry"
        compo.reduced formula (str): The formula of the identified
                                     network
        list(og_cluster) (list): The list of sites associated with
                                 the network, relative to the original
                                 structure object
    """
    
    
    # Entry is a list of the necessary components of the TSA. 
    # Makes it easier to parallelize
    structure, tol, mp_id = entry
    norm_tol=tol-1

    s = copy.deepcopy(structure)
    heterogeneous = False
    heterogeneousSS = False

    # Distance matrix (rowA, columnB) shows distance between
    # atoms A and B, taking PBCs into account.
    
    binary_matrix = getDistMat(s,norm_tol)
    cluster = buildNetwork(binary_matrix,seed_index)

    compo = Composition.from_dict(Counter([s[l].specie.name for 
                                           l in list(cluster)]))
    if compo.reduced_formula != s.composition.reduced_formula:
        # i.e. the cluster does not have the same composition
        # as the overall crystal; therefore there are other
        # clusters of varying composition.
        heterogeneous = True
    og_cluster = set()
    og_cluster.update(cluster)
    old_cluster_size = len(cluster)
    # Increase structure to determine dimensionality


    s.make_supercell(supercell)
    seed_index*=supercell**3



    binary_matrix = getDistMat(s,norm_tol)
    cluster = buildNetwork(binary_matrix,seed_index)


    if cluster!=set():
        new_cluster_size = len(cluster)
        # Get ratio of original and final cluster lengths
        scale = new_cluster_size/old_cluster_size
        compo = Composition.from_dict(Counter([s[l].specie.name for l in
                                    list(cluster)]))
        if compo.reduced_formula != s.composition.reduced_formula:
            # i.e. the cluster does not have the same composition
            # as the overall crystal; therefore there are other
            # clusters of varying composition.
            heterogeneousSS = True
        motiif = getDim(scale,supercell)
    if heterogeneous or heterogeneousSS:
        motiif += "_heter"

    if returnSS:
        return [motiif,mp_id,tol,compo.reduced_formula,list(cluster)]
    else:
        return [motiif,mp_id,tol,compo.reduced_formula,list(og_cluster)]



def checkHeter(structure,batch,supercell=2):
    """
    A function to search for heterostructures. If the "getStructureType"
    does not report the network as a heterostructure, the input network
    will be returned. If a heterostructure is identified, then the TSA
    will be performed on the structure until all atoms are associated 
    with an atomic network. Any dupicate structures will need to be 
    removed through post-processing
    
    inputs
    --------
        structure (Structure): pymatgen Structure object of the original 
                               crystal
        batch (list)         : the output of "getStructureType"

    returns
    --------
        batches (list) : list of "getStructureType" outputs for each 
                         network identified
        uniqueStructs (list) : A list of Structure objects that represent
                               each network found in the system
    """
    s = copy.deepcopy(structure)
    # Check if structure is heterostructure
    if 'heter' in batch[0]:
        i=0
        found = list(batch[4])
        all_structs=[s.from_sites([s.sites[n] for n in batch[4]])]
        batches = [batch]
        # Loop through all sites
        for i in range(s.num_sites):
            # Skip sites that have already been accounted for
            if i not in found: 
                s = copy.deepcopy(structure)
                nBatch = getStructureType([s,batch[2],batch[1]],
                                          seed_index=i,supercell=supercell)
                if nBatch[4][0] not in found:
                    found.extend(nBatch[4])
                    all_structs.append(structure.from_sites([s.sites[n] for 
                                                             n in nBatch[4]]))
                batches.append(nBatch)
        return(batches,all_structs)
    else:
        return([batch],[s.from_sites([s.sites[n] for n in batch[4]])])


def getSpecDict(structure):
    # Get dictionary of species in s [Structure object]

    spec_dict={}
    for site in structure.sites:
        if site.specie not in spec_dict:
            spec_dict[site.specie]=1
        else:
            spec_dict[site.specie]+=1


    # Create a starting point for the looping
    refList = [site.specie,spec_dict[site.specie]]
    for key in spec_dict.keys():
        if spec_dict[key]<refList[1]:
            refList = [key,spec_dict[key]]

    # Get the most infrequent specie
    ref_spec=refList[0]


    return(ref_spec)


def magni(vector):
    """
    Helper funciton for getting the magnitude of a vector
    
    inputs
    --------
        vector (array 1xN or Nx1):  An array of numbers representing 
                                    a vetor in space

    returns
    --------
        float representing the magnitude of the input vector
        
    """
    return(np.linalg.norm(vector))


def cleaveSurfAtom(entry,max_bonds=1,supercell=2,group_structs=True):
    """
    An algorithm to cleave a surface from a fully periodic crystal. This 
    iteration uses the "periodic atom" approach, as described in Paul et al.
    in 2020. If this algorithm is used, please cite:
        
        https://arxiv.org/abs/2002.00903
    
    inputs
    --------
        entry (list): A set of components necessary for the TSA.
                      Makes it easier to parallelize with this as
                      the input
                      --structure (Structure): pymatgen Structure object
                      --tol (float): The scaling for the atomic bonds
                      --mp_id (str): The label for the entry, commonly
                                     the MaterialsProject ID
        max_bonds (int):    The maximum number of bonds to cleave with this
                            algorithm. 
                            
        supercell (int):    The size supercell to use for building the 
                            bonded atomic network. Typically does not 
                            matter, except in the case of a bipartide
                            crystal. In such cases, even or odd sized
                            supercells will result in the TSA identifying 
                            a clear dimensionality or "mixed" dimensionality,
                            not respectivel.
        group_structs (bool):   Whether to group the surfaces cleaved based 
                                on symmetry, or give the full list of surfaces
                                cleaved. The more symmetric the crystal, the
                                more duplicate surfaces will be present                                    

    returns
    --------
        cleaved surfaces (list): The list of structure objects representing
                                 the cleaved structures. If it is empty,
                                 no surfaces could be cleaved
        
    """
    
    
    struct = copy.deepcopy(entry[0])
    results = getStructureType(entry,supercell=supercell,returnSS=True)
    if results[0]=='conventional':
        struct = copy.deepcopy(entry[0])
        og_binary_matrix = getDistMat(struct,entry[1]-1)
        og_num_bonds = sum(sum(np.array(og_binary_matrix)))
        struct.make_supercell(supercell)
        binary_matrix= getDistMat(struct,entry[1]-1)
        bonds = []
        for i in range(len(og_binary_matrix)):
            for pair in [(i,j) for j in range(i+1,len(og_binary_matrix)) 
                         if og_binary_matrix[i][j]==1]:
                bonds.append(pair)
        allCombos = []
        combNum = 0
        for i in range(max_bonds+1):
            for com in list(itertools.combinations(bonds,i)):
                allCombos.append(com)
                combNum+=1

        combos = allCombos
        jjj=0
        all_structs = []
        for combo in combos:
            broken=0
            jjj+=1
            modified_matrix = np.array(binary_matrix)
            for pair in combo:
                i,j = pair
                i=i*supercell**3
                j=j*supercell**3
                for shift in range(supercell**3):
                    for shift2 in range(supercell**3):
                        modified_matrix[i+shift][j+shift2]=0
                        modified_matrix[j+shift][i+shift2]=0

            new_num_bonds=sum(sum(modified_matrix))
            broken=int(og_num_bonds-new_num_bonds)  
            seed_index=0
            old_cluster_size=len(buildNetwork(binary_matrix,seed_index))/supercell**3
            cluster = buildNetwork(modified_matrix,seed_index)
            hetero=False
            if cluster!=set():
                scale = len(cluster)/old_cluster_size
                compo = Composition.from_dict(Counter([struct[l].specie.name 
                                                 for l in list(cluster)]))
                if compo.reduced_formula != struct.composition.reduced_formula:
                    # i.e. the cluster does not have the same composition
                    # as the overall crystal; therefore there are other
                    # clusters of varying composition.
                    hetero = True
                motiif = getDim(scale,supercell)

            if not hetero:
                if motiif=='layered':
                    cluster_sites = [struct.sites[n] for n in cluster]
                    all_structs.append(struct.from_sites(cluster_sites))

        if group_structs:
            matched = [x[0] for x in 
                       StructureMatcher(stol=1E-6,primitive_cell=False,
                       scale=False).group_structures(all_structs)]
        else:
            matched=all_structs
        return(matched)    

    else:
        print('Material is does not have a 3D motiif')
        print('Try increasing radii tolerance if appropriate')
        return([])

def getBondVectors(struct,tol,prec):
    
    """
    Helper function for "cleaveSurfBond" for getting the list of 
    bonds which are parallel and of identical magnitude
    
    inputs
    --------
        struct (Structure): pymatgen Structure object
        
        tol (float): The scaling for the atomic radii when building
                     the atomic networks

        prec (float):   The precision to compare magnitude of vectors
                        representing the bonds in the system                                      

    returns
    --------
        bond_dir (dict): Dictionary containing list of all bonds in the 
                         structure
        
    """    
    
    
    binary_matrix= getDistMat(struct,tol)
    bond_dir = {}
    distance_matrix = struct.distance_matrix
    lattice = np.array(struct.lattice.as_dict()['matrix'])
    iterations = list(itertools.product([1,0,-1],repeat=3))
    # Loop over list of atoms
    for i in range(len(binary_matrix)):
        for j in range(i+1,len(binary_matrix)):
            # Proceed if the entries are listed as "bonded"      
            if binary_matrix[i][j]==1:                 
                s1 = struct.species[i]
                s2 = struct.species[j]
                # Organize dictionary so it is always in order of increasing
                # atomic number
                if s1.number>s2.number:
                    s1 = struct.species[j]
                    s2 = struct.species[i] 
                if s1 not in bond_dir:
                    bond_dir[s1]={}
                if s2 not in bond_dir[s1]:
                    bond_dir[s1][s2]=[]
                valid_vs = []
                
                # Get the vector between atomic positions
                
                bond_vector = np.array(struct.sites[j].coords-
                                          struct.sites[i].coords) 
                
                # The positions of the atoms may not be in the right locations
                # to be the minimum distance from each other. As a result,
                # a translation is applied to the resulting "bond vector" 
                # (alternatively, one of the atoms is translated)
                for shift in iterations:
                    bondShift = bond_vector + np.dot(lattice.T,shift)
                    if abs(distance_matrix[i][j]-magni(bondShift))<=prec:
                            valid_vs.append(bondShift)
                            break
                # See if the vector is already present in the collection of 
                # vectors. If so, add the coordinates to the entry. Else,
                # create a new entry for the direction of the bond.
                for v in valid_vs:
                    if np.any([magni(v-x[0])<=prec for x in bond_dir[s1][s2]]):
                        for k in range(len(bond_dir[s1][s2])):
                            if magni(v-bond_dir[s1][s2][k][0])<=prec:
                                bond_dir[s1][s2][k][1].append([i,j])
                                break
                    
                    else:
                        bond_dir[s1][s2].append([v,[[i,j]]])
    return(bond_dir)

def cleaveSurfBond(entry,max_bonds=1,supercell=2,group_structs=True,prec=1E-4):
    
    """
    An algorithm to cleave a surface from a fully periodic crystal. This 
    iteration uses the "periodic bond" approach, as described in Paul et al.
    in 2020. If this algorithm is used, please cite:
        
        https://arxiv.org/abs/2002.00903
    
    inputs
    --------
        entry (list): A set of components necessary for the TSA.
                      Makes it easier to parallelize with this as
                      the input
                      --structure (Structure): pymatgen Structure object
                      --tol (float): The scaling for the atomic bonds
                      --mp_id (str): The label for the entry, commonly
                                     the MaterialsProject ID
                                     
        max_bonds (int):    The maximum number of bonds to cleave with this
                            algorithm. 
                            
        supercell (int):    The size supercell to use for building the 
                            bonded atomic network. Typically does not 
                            matter, except in the case of a bipartide
                            crystal. In such cases, even or odd sized
                            supercells will result in the TSA identifying 
                            a clear dimensionality or "mixed" dimensionality,
                            not respectively.
                            
        group_structs (bool):   Whether to group the surfaces cleaved based 
                                on symmetry, or give the full list of surfaces
                                cleaved. The more symmetric the crystal, the
                                more duplicate surfaces will be present
                                
        prec (float):   The precision to compare magnitude of vectors
                        representing the bonds in the system                                      

    returns
    --------
        cleaved surfaces (list): The list of structure objects representing
                                 the cleaved structures. If it is empty,
                                 no surfaces could be cleaved
        
    """
    
    
    struct = copy.deepcopy(entry[0])
    results = getStructureType(entry,supercell=supercell,returnSS=True)
 
    # Proceed only if the structure is classified as periodic
    # in all directions
    if results[0]=='conventional':
        struct.make_supercell(supercell)
        binary_matrix= getDistMat(struct,entry[1]-1)
        og_num_bonds = sum(sum(np.array(binary_matrix)))/2
        
        # Get dictionary of directional bonds in the system,  
        # and the associated atomic species
        bond_dir = getBondVectors(struct,entry[1]-1,prec)

                        
        # Create the list of bonds to be broken
        all_structs=[]
        combos=[]
        for s1 in bond_dir:
            for s2 in bond_dir[s1]:
                for cleave in bond_dir[s1][s2]:    
                        combos.append(cleave[1])
                        
        # Create pairings of bonds to be broken, up to 
        # max_bonds number of bonds
                        
        final_combos=[]
        for i in range(1,max_bonds+1):
            for mix in list(itertools.combinations(combos,max_bonds)):
                    final_combos.append(mix)
        seed_index=0
        old_cluster_size=len(buildNetwork(binary_matrix,seed_index))/supercell**3
        for combo in final_combos:
            modified_matrix = np.array(binary_matrix)
            for sett in combo:
                for pair in sett:
                    i,j = pair
                    modified_matrix[i][j]=0
                    modified_matrix[j][i]=0
            new_num_bonds=sum(sum(modified_matrix))/2
            
            # Number of bonds broken in the search. Not necessarily
            # the number of bonds broken to cleave the surface
            
            broken=int(og_num_bonds-new_num_bonds)
            
            cluster = buildNetwork(modified_matrix,seed_index)
            hetero=False
            if cluster!=set():
                scale = len(cluster)/old_cluster_size
                compo = Composition.from_dict(Counter([struct[l].specie.name 
                                             for l in list(cluster)]))
                if compo.reduced_formula != struct.composition.reduced_formula:
                    # i.e. the cluster does not have the same composition
                    # as the overall crystal; therefore there are other
                    # clusters of varying composition.
                    hetero = True
                motiif = getDim(scale,supercell)

            if not hetero:
                if motiif=='layered':
                    cluster_sites = [struct.sites[n] for n in cluster]
                    all_structs.append(struct.from_sites(cluster_sites))

        if group_structs:
            matched = [x[0] for x in 
                       StructureMatcher(stol=1E-6,primitive_cell=False,
                       scale=False).group_structures(all_structs)]
        else:
            matched=all_structs
        return(matched)    


    else:
        print('Material is does not have a 3D motiif')
        print('Try increasing radii tolerance if appropriate')
        return([])
    
    
    
def getAtomImages(struct,ref_site,supercell=2,prec=1E-4):
    
    """
    
    Helper funciton for getting the periodic 
    images of an atom in the supercell
    
    inputs
    --------
        struct (Structure): pymatgen Structure object
        
        ref_sites (list, 3x1):  the cartesian coordinates of the reference
                                species
        
        supercell (int):  The size of the supercell to generate

        prec (float):   The precision to compare magnitude of vectors
                        representing the bonds in the system 

    returns
    --------
        list of indexes in "struct" which are periodic images of 
        ref_site
    
    
    """
    
    ref_latt=np.array(struct.lattice.as_dict()['matrix'])/supercell
    lattice_shifts = list(itertools.product([1,0,-1],repeat=3))
    # Get coordinates of periodic images of the reference atom
    ref_coords = [np.array(ref_site+np.dot(ref_latt.T,x)) for 
                  x in lattice_shifts]
    # Get index of periodic indexes of the reference atom
    periodic_sites=[atom_ind for atom_ind in range(struct.num_sites) if \
                    np.any([np.linalg.norm(struct.sites[atom_ind].coords-\
                    ref_site)<=prec for ref_site in ref_coords])]
    return(periodic_sites)





def getNewLattice(entry,dim,prec=1E-4,seed_index=0,supercell=2,c_mag=50):
    
    """
    Helper function for "cleaveSurfBond" for getting the list of 
    bonds which are parallel and of identical magnitude
    
    inputs
    --------
       entry (list): A set of components necessary for the TSA.
                      Makes it easier to parallelize with this as
                      the input
                      --structure (Structure): pymatgen Structure object
                      --tol (float): The scaling for the atomic bonds
                      --mp_id (str): The label for the entry, commonly
                                     the MaterialsProject ID
        
        dim (int):    Number of periodic directions in structure
        
        prec (float):       The precision to compare magnitude of vectors
                            representing the bonds in the system     
                        
        seed_index (int):   The site to use as the starting point for the
                            TSA. Typically does not impact the results, but
                            will if the structure is a bipartide or has
                            mixed dimensionality    
                            
        supercell (int):    The supercell size to generate for 
                            periodic networks
                            
        c_mag (float):      The magnitude to make the non-periodic vectors
        
        
    returns
    --------
        latt_attempt (list, 3X3): New lattice for "structure"
        
    """   
    
    
    
    # Get network of atoms within the supercell

    # Reset structure object, as "getStructureType" turned it into a supercell
    structure = entry[0]
    # Shift atoms so that the refernce atom is at the center of the cell
    structure.translate_sites(indices=range(structure.num_sites),
                              vector=-1*structure.sites[0].frac_coords+
                              [.02,.02,.02],frac_coords=True)
 
    og_structure = copy.deepcopy(structure)
    
    # Replace structure used for TSA with centered structure
    entry[0]=structure
    batch = getStructureType(entry,seed_index=seed_index,
                             supercell=supercell,returnSS=True)

    # Get position of reference atom, and periodic images
    ref_site = np.array(structure.sites[0].coords)


    # Get structure object containing a single atomic network
    structure.make_supercell(supercell)
    cluster_sites = [structure.sites[n] for n in batch[4]]
    new_struct = structure.from_sites(cluster_sites)

    # Get the list of reference atoms which are in the
    # isolated atomic network
    in_network = [new_struct.sites[x].coords for x in getAtomImages(new_struct,ref_site)]

    # Ensure all fractional coordiantes are <1

    for i in range(og_structure.num_sites):
        if np.any([abs(x)>=1 for x in og_structure.sites[i].frac_coords]):
            tVector = [x-x%1 for x in og_structure.sites[i].frac_coords]
    #        print(ogStructure.sites[i].frac_coords)
    #        print(tVector)
            og_structure.translate_sites(indices=i,vector=-1*tVector)



    # Get the new lattice matrix
    latt_attempt = genLattice(og_structure,in_network,dim)
    
    # If a valid lattice matrix isn't found, it is likely due to an issue
    # with original lattice vectors not being orthogonal. 
    # Shift atomic positions a bit until this problem doesn't occur
    o_shift=0.05
    o_shift_count=float(o_shift)
    while latt_attempt==[] and o_shift_count<1:
        print(o_shift_count)
        og_structure.translate_sites(indices=range(og_structure.num_sites),
                                    vector=[o_shift,o_shift,o_shift],
                                    frac_coords=True)
        og_structure.to('POSCAR','POS_'+str(round(o_shift_count,2)))
        o_shift_count+=.05
        latt_attempt = genLattice(og_structure,in_network,dim)
    return(latt_attempt)


def calcCrossMag(v1,v2):
    """
    Helper function for getting the magnitude of a cross product.
    Used primarily to determine if two vectors are parallel
    
    inputs
    --------
        v1,v2 (list,list): Two lists representing vectors
        
        
    returns
    --------
        the magnitude of the cross product of v1 and v2
    
    """
    # Calculate the magnitude of cross product of two vectors

    return(abs(np.linalg.norm(np.cross(v1,v2))))


def getTranslation(fracs):
    """
    Get a translation vector such at all atoms are as close
    to the origin as possible. Used in the "genLattice" function 
    to check if the new unit cell captures all the atoms within it.
    
    inputs
    --------
        fracs (list):  Nx3 matrix of the fractional coordinates of 
                       a crystal. The first, second, and third vector
                       should be the x, y, and z fractional coordinates
                       in the crystal
    returns
    --------
        shift_vector (list) :   vector that should be subtracted from all
                                fractional coordinates of each site
    
    """
    
    
    
    # Determine whether the shift needs to be from inf to 0 
    # or from -inf to 0
    
    # Along all x fractionals
    if abs(max(fracs[0]))>=abs(min(fracs[0])):
        minX = min([x for x in fracs[0] if x>0])
    else:
        minX = min([x for x in fracs[0] if x<0])
        
    # Along all y fractionals
    if abs(max(fracs[1]))>=abs(min(fracs[1])):
        minY = min([x for x in fracs[1] if x>0])
    else:
        minY = min([x for x in fracs[1] if x<0])
        
    # Along all z fractionals
    # Need to consider all atoms lying in a single
    # plane (e.g. graphene), thus the final "else"
    # statement
    if abs(max(fracs[2]))>abs(min(fracs[2])):
        minZ = min([x for x in fracs[2] if x>0])
    elif abs(max(fracs[2]))<abs(min(fracs[2])):
        minZ = min([x for x in fracs[2] if x<0])
    else:
        minZ = max(fracs[2])

    shift_vector = np.array([minX,minY,minZ])
   
    return(shift_vector)
    

def genLattice(structure,in_network,dimension,prec=1E-4,
               seed_index=0,supercell=2,c_mag=50):
    
    """
    Generate a new lattice for the low-dimensional material. Creates a lattice
    such that the vectors representing non-periodic directions are orthogonal 
    to the periodic directions
    
    inputs
    --------
        structure (Structure): pymatgen Structure object
        
        in_network (list):  List of atoms which are periodic images of 
                            each other
        
        dimension (int):    Number of periodic directions in structure
        
        prec (float):       The precision to compare magnitude of vectors
                            representing the bonds in the system     
                        
        seed_index (int):   The site to use as the starting point for the
                            TSA. Typically does not impact the results, but
                            will if the structure is a bipartide or has
                            mixed dimensionality    
                            
        supercell (int):    The supercell size to generate for 
                            periodic networks
                            
        c_mag (float):      The magnitude to make the non-periodic vectors
            
    returns
    --------
        new lattice (list): List of lists representing the new lattice for
                            the low-dimensional material
        
    """

    # Generate vectors in plane/line, relative to
    # the first atom in the network of atoms
    vec_list = {}
    for i in range(len(in_network)-1):
        vec_list[i]=[]
        for j in range(i,len(in_network)):
            if i!=j:
                vec_list[i].append(in_network[j]-in_network[i])
    print(vec_list)
    vect_and_magnitude = []
    if dimension==2:
        for i in vec_list:
            for sett in [[[x,y],magni(np.cross(x,y))] for [x,y] in list(itertools.combinations(vec_list[i],dimension))]:
                vect_and_magnitude.append(sett)
    elif dimension==1:
        for i in vec_list:
            for sett in [[x,magni(x)] for x in list(itertools.combinations(vec_list[i],dimension))]:
                vect_and_magnitude.append(sett)
        
    vectors = [np.array(pair) for [pair,magnitude] in sorted(vect_and_magnitude,key=operator.itemgetter(1))]
    coords = [x.coords for x in structure.sites]
    return_structure=False
    
    # Attempt all mixtures of vectors as
    # potential lattice vectors
    for lat_vectors in vectors:
        # Create lattice matrix to fit atomic coordinates against
        
        # In 2D
        if dimension==2:
                scaleC=c_mag/magni(lat_vectors[0])/magni(lat_vectors[1])
                latt_attempt = np.array([lat_vectors[0],lat_vectors[1],\
                    np.cross(lat_vectors[0],lat_vectors[1])*scaleC])
                
        # In 1D
        elif dimension==1:
            unitV = lat_vectors[0]/np.linalg.norm(lat_vectors[0])
            if unitV[0]==0:
                perp1 = [1,0,0]
            elif unitV[1]==0:
                perp1 = [0,1,0]
            elif unitV[2]==0:
                perp1 = [0,0,1]
            else:
                perp1 = [1,1,-1*(unitV[0]+unitV[1])/unitV[2]]
            perp1 = perp1/np.linalg.norm(perp1)*c_mag
            perp2 = np.cross(unitV,perp1)
            perp2 = perp2/np.linalg.norm(perp2)*c_mag
            latt_attempt   = np.array([lat_vectors[0],perp1,perp2])
        # Fit atomic sites to new lattice
        temp_fracs = np.linalg.solve(latt_attempt.T,np.array(coords).T)

        # Make list of all fractional positions, ignoring
        # which axis
        new_fracs = list([list(x) for x in temp_fracs.T])
        temp_cat = new_fracs[0]+new_fracs[1]+new_fracs[2]
        
        # Ensure cartesian coordinate of all 
        # atoms can be described using no more 
        # one of each unit vector
        if not np.any([abs(round(x,2))>1 for x in temp_cat]):
            return_structure=True
            break

    # If match found
    if return_structure:
        return(latt_attempt)
    # If no match found
    else:
        return([])


def alignMono(entry,prec=1E-4,seed_index=0,supercell=2,
              c_mag=50,dist_from_plane=0):
    
    """
    Align a 2D material such that the 'c' vector is perpendicular to the
    in-plane lattice vectors
    
    inputs
    --------
       entry (list): A set of components necessary for the TSA.
                      Makes it easier to parallelize with this as
                      the input
                      --structure (Structure): pymatgen Structure object
                      --tol (float): The scaling for the atomic bonds
                      --mp_id (str): The label for the entry, commonly
                                     the MaterialsProject ID

        
        prec (float):       The precision to compare magnitude of vectors
                            representing the bonds in the system     
                        
        seed_index (int):   The site to use as the starting point for the
                            TSA. Typically does not impact the results, but
                            will if the structure is a bipartide or has
                            mixed dimensionality    
                            
        supercell (int):    The supercell size to generate for 
                            periodic networks
                            
        c_mag (float):      The magnitude to make the non-periodic vectors
        
        dist_from_plane (float): Maximum distance an atom can be from the
                                 plane parallel to the monolayer. Is relevant 
                                 when the atoms in the monolayer are spread 
                                 across periodic boundary conditions in the 
                                 unit cell
                                 
        
        
    returns
    --------
        list1 (list): -fractional coordinates in new lattice
                     -new lattice (a,b,c)
                     -species associated with each site
        list1 (list): -fractional coordinates in new lattice
                     -new lattice (b,a,c)
                     -species associated with each site
    
    """


    # Keep original copy of structure
    s = copy.deepcopy(entry[0])


    new_latt = getNewLattice(entry,dim=2,prec=prec,seed_index=seed_index,
                          supercell=supercell,c_mag=c_mag)


    # Identify plane to translate atoms towards
    plane = Plane(Point3D(s.sites[seed_index].coords),
                  normal_vector=new_latt[2])
    
    # Create list of translationss
    trans = list(itertools.product([1,-1,0],repeat=3))

    lat = s.lattice.as_dict()['matrix']
    final_sites = []
    i=0
    
    # Ensure that the atoms are nearby each other
    for site in s.sites:
        point = Point3D(site.coords)
        if plane.distance(point)<dist_from_plane:
            final_sites.append(site.coords)
        else:
            news = []
            
            # translate atomic sites to see which position is closest to plane
            for t in trans:
                point = Point3D(site.coords+np.dot(np.transpose(lat),t))
                news.append([float(plane.distance(point)),t])
            news.sort(key = lambda x:x[0])
            final_sites.append(site.coords+
                               np.dot(np.transpose(lat),news[0][1]))
        i+=1
        
    # Create new lattice matricies
    lat1 = np.array([new_latt[0],new_latt[1],new_latt[2]])
    lat2 = np.array([new_latt[1],new_latt[0],new_latt[2]])

    # Generate atomic fractions
    new_fracs1 = np.linalg.solve(lat1.T,np.array(final_sites).T).T
    new_fracs2 = np.linalg.solve(lat2.T,np.array(final_sites).T).T
    species = s.species
    return([species,new_fracs1,lat1],[species,new_fracs2,lat2])

def alignChain(entry,prec=1E-4,seed_index=0,supercell=2,
               c_mag=50,dist_from_line=0):
        
    """
    Align a 2D material such that the 'c' vector is perpendicular to the
    in-plane lattice vectors
    
    inputs
    --------
       entry (list): A set of components necessary for the TSA.
                      Makes it easier to parallelize with this as
                      the input
                      --structure (Structure): pymatgen Structure object
                      --tol (float): The scaling for the atomic bonds
                      --mp_id (str): The label for the entry, commonly
                                     the MaterialsProject ID

        
        prec (float):       The precision to compare magnitude of vectors
                            representing the bonds in the system     
                        
        seed_index (int):   The site to use as the starting point for the
                            TSA. Typically does not impact the results, but
                            will if the structure is a bipartide or has
                            mixed dimensionality    
                            
        supercell (int):    The supercell size to generate for 
                            periodic networks
                            
        c_mag (float):      The magnitude to make the non-periodic vectors
        
        dist_from_line (float):  Maximum distance an atom can be from the
                                 line parallel to the nanowire. Is relevant 
                                 when the atoms in the nanowire are spread 
                                 across periodic boundary conditions in the 
                                 unit cell
                                 
        
        
    returns
    --------
        list (list): -fractional coordinates in new lattice
                     -new lattice (a,b,c)
                     -species associated with each site
    
    """

    new_struct = copy.deepcopy(entry[0])

    new_latt = getNewLattice(entry,1,prec,seed_index,supercell,c_mag)
    print(new_latt)
    v1,v2,perp=new_latt
    new_latt=np.array(new_latt)
    line1 = Line3D(new_latt[0],[0,0,0])
    line = line1.parallel_line(new_struct.sites[0].coords)
    trans = list(itertools.product([1,-1,0],repeat=3))
    print('WORKING1')
    lat = np.array(new_struct.lattice.as_dict()['matrix'])
    final_sites = []
    i=0
    i=-1
    for site in [x.coords for x in new_struct.sites]:
        print(i,new_struct.num_sites)
        i+=1
        point = Point3D(site)
        if line.distance(point)<dist_from_line:
            final_sites.append(site)
        else:
            news = []
            for t in trans:
                point = Point3D(site+np.dot(lat.T,t))
                news.append([float(line.distance(point)),t])
            news.sort(key = lambda x: x[0])
            final_sites.append(site+np.dot(lat.T,news[0][1]))


    #new_latt = np.array([,perp2])
    new_fracs = np.linalg.solve(new_latt.T,np.array(final_sites).T).T
    species = new_struct.species
    return([species,new_fracs,new_latt])


def getAngle(v1,v2,prec=1E-6):

    """
    Helper function to get angle between two vectors
    
    inputs
    --------
        v1,v2 (list,list): two vectors
        
    outputs
    --------
        angle (float, radians): angle between v1,v2
    
    """
    
    return(math.acos((np.dot(v1,v2))/np.linalg.norm(v1)/np.linalg.norm(v2)))


def makeNewPos(specs,frac_coords,new_latt,dim):
    
    
    """
    Make a new POSCAR file using input species, coordinates, and lattice.
    Used in conjunction with functions "getNewLattice","alignMono",
    "alignChains"
    
    inputs
    --------
        specs (list):       List of species associated with each coordinate
        
        frac_coords (list): List of fractional coordinates of each species,
                            relative to the new lattice
                            
        new_latt (list):    List of vectors representing the new lattice 
                            vectors of the low-dimensional material
                            
        dim (int):  Number of periodic directions in the crystal
        
    outputs
    --------
        new_struct (Structure): Structure object for low-dimensional material
                                with non-periodic directions being 
                                orthogonal to the periodic directions
    
    """

    
    a,b,c = magni(new_latt[0]),magni(new_latt[1]),magni(new_latt[2])
    if dim==2:
        ang = getAngle(new_latt[0],new_latt[1])/2
        new_latt = [[np.cos(ang)*a,-np.sin(ang)*a,0],
                    [np.cos(ang)*b,np.sin(ang)*b,0],[0,0,c]]
    elif dim==1:    
        new_latt = [[a,0,0],[0,b,0],[0,0,c]]
    i=0
    new_sites = []
    for site in frac_coords:
        p = PeriodicSite(species=Element(specs[i]),
                         lattice = Lattice(new_latt),
                         coords=np.dot(site,new_latt),
                         coords_are_cartesian=True)
        new_sites.append(p)
        i+=1

    new_struct = Structure.from_sites(new_sites)
    return(new_struct)


def reduceScale(structure,scale,dim):
    
    """
    Attempt to reduce the structure from a supercell
    
    inputs
    --------
        structure (Structure):  Pymatgen structure object to reduce
        
        scale (float): How much to scale the periodic lattice vectors by
        
        dim (int):  The number of periodic directions in "structure"
        
        
    outputs
    --------
        new_struct (Structure): Structure object for low-dimensional material
                                with non-periodic directions being 
                                orthogonal to the periodic directions
    
    """

    structure.translate_sites(indices=range(structure.num_sites),
                            vector=-1*structure.sites[0].frac_coords+[0,0,.5])

    specs = structure.species
    cart = [x.coords for x in structure.sites]

    lat = np.array(structure.lattice.as_dict()['matrix'])
    lat[0]*=scale
    if dim==2 or dim==3:
        lat[1]*=scale
    if dim==3:
        lat[2]*=scale
    fracs = np.array(np.linalg.solve(lat.T,np.array(cart).T).T)
    specs   = []
    u_cart   = []
    i=0
    for frac in fracs:
        if frac[0]<1 and frac[1]<1 and frac[2]<1:
                specs.append(structure.species[i])
                u_cart.append(cart[i])
        i+=1
    
    new_sites = []
    i=0
    for site in u_cart:
        p = PeriodicSite(atoms_n_occu=Element(specs[i]),
                         lattice = Lattice(lat),coords=site,
                         coords_are_cartesian=True)
        new_sites.append(p)
        i+=1


    new_struct = Structure.from_sites(new_sites)

    return(new_struct)


def getPlane(entry): 
     
    """
    Get the (X,Y,Z) of the plane of the monolayer in the original lattice
    
    inputs
    --------
       entry (list): A set of components necessary for the TSA.
                      Makes it easier to parallelize with this as
                      the input
                      --structure (Structure): pymatgen Structure object
                      --tol (float): The scaling for the atomic bonds
                      --mp_id (str): The label for the entry, commonly
                                     the MaterialsProject ID
    outputs
    --------
        fracs (list): List of x,y,z intercepts of the monolayer
    
    """

    
    
    a,b,c = getNewLattice(entry,2)
    a_vector = np.linalg.solve(np.array(entry[0].lattice.as_dict()['matrix']).T,a)
    b_vector = np.linalg.solve(np.array(entry[0].lattice.as_dict()['matrix']).T,b)
    fracs = np.cross(a_vector,b_vector)
    fracs /= min([x for x in fracs if abs(x)>1E-4])
    
    return(fracs)

