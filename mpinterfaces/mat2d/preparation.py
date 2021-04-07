from __future__ import print_function, division, unicode_literals

import copy
import numpy as np
import math
import itertools


from collections import Counter
from sympy import Point3D, Plane,Line3D


from operator import itemgetter
from pymatgen import Structure, Element, Composition
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import _pt_data

__author__ = "Joshua Paul"
__copyright__ = "Copyright 2021, Henniggroup"
__maintainer__ = "Joshua Paul"
__email__ = "joshua.thomas.paul@gmail.com"
__status__ = "Production"
__date__ = "March 2, 2021"


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
        print('WARNING: Increase in radius is greater than 100%')
        
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

def getStructureType(entry, seed_index=0, supercell = 2,return_SS=False):
    
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
    heterogeneous_SS = False

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

    s = copy.deepcopy(structure)
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
            heterogeneous_SS = True
        motiif = getDim(scale,supercell)
    if heterogeneous or heterogeneous_SS:
        motiif += "_heter"

    if return_SS:
        return [motiif,mp_id,tol,compo.reduced_formula,list(cluster)]
    else:
        return [motiif,mp_id,tol,compo.reduced_formula,list(og_cluster)]


def checkHeter(structure,batch,supercell=2,try_all=False,return_SS=False):
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
        supercell (int)      : how large of a supercell to use when
                               searching for networks
        try_all (Boolean)    : whether to search across all atoms,
                               regardless of dimensionality of input
                               structure
        return_SS (Boolean)  : whether to return the atomic network from the
                               supercell instead of primitive cell
                           
    returns
    --------
        batches (list) : list of "getStructureType" outputs for each 
                         network identified
        uniqueStructs (list) : A list of Structure objects that represent
                               each network found in the system
    """
    s = copy.deepcopy(structure)
    # Check if structure is heterostructure
    if 'heter' in batch[0] or try_all:
        i=0
        
        found=[]
        all_structs=[]
        batches=[]
        # Loop through all sites
        for i in range(s.num_sites):
            
            # Skip sites that have already been accounted for
            if i not in found: 

                s = copy.deepcopy(structure)
                nBatch = getStructureType([s,batch[2],batch[1]],
                                          seed_index=i,supercell=supercell,
                                          return_SS=return_SS)

                if np.any(x not in found for x in nBatch[4]):
                    found.extend(nBatch[4])
                    s = copy.deepcopy(structure)
                    s.make_supercell(supercell)
                    all_structs.append(Structure.from_sites([s.sites[n] for 
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
    
    # If the crystal is 3D
    if results[0]=='conventional':
        struct = copy.deepcopy(entry[0])
        og_binary_matrix = getDistMat(struct,entry[1]-1)
        og_num_bonds = sum(sum(np.array(og_binary_matrix)))
        struct.make_supercell(supercell)
        binary_matrix= getDistMat(struct,entry[1]-1)
        bonds = []
        
        # Get list of bonded atoms
        for i in range(len(og_binary_matrix)):
            for pair in [(i,j) for j in range(i+1,len(og_binary_matrix)) 
                         if og_binary_matrix[i][j]==1]:
                bonds.append(pair)
        allCombos = []
        combNum = 0
        
        # Get list of all combinations of bonds
        for i in range(max_bonds+1):
            for com in list(itertools.combinations(bonds,i)):
                allCombos.append(com)
                combNum+=1

        combos = allCombos
        jjj=0
        all_structs = []
        
        # For each bond combination
        for combo in combos:
            broken=0
            jjj+=1
            modified_matrix = np.array(binary_matrix)
            for pair in combo:
                i,j = pair
                i=i*supercell**3
                j=j*supercell**3
                # Break bonds in the loop
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
            # If the new set of atoms is not empty
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

            # If formula of new network matches the original cell
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
    
    Helper function for getting the periodic 
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



def getUniqueCount(specs):
    '''
    
    Helper function for counting the number of instances
    of a species
    
    inputs
    ---------
        specs (list): List of pymatgen element objects
        
        
    returns
    ---------
        a dictionary with keys being pymatgen element objects
        and the output being the number of instances of that element
    
    '''
    
    d = {}
    for entry in specs:
        if entry not in d:
            d[entry]=1
        else:
            d[entry]+=1
            
    return(d)



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
    

    
def reduceCoords(structure):
    
    '''
    helper function for reducing fractional coordinates
    
    
    '''
    
    # Ensure all fractional coordiantes are <1
    sites = structure.sites
    for i in range(structure.num_sites):
        if np.any([round(abs(x),4)>=1 for x in sites[i].frac_coords]):
            tVector = np.array([x-round(x,4)%1 for x in sites[i].frac_coords])


            structure.translate_sites(indices=i,vector=-1*tVector,
                                        frac_coords=True)


    return(structure)
    

def getVectorsRough(structure, supercell, prec,dimension,seed_index=0):
    '''
    
    Helper function to get the in-plane lattice vectors 
    in a more complete, but significantly less efficient way.
    Typically not needed, except in extreme cases.
    
        inputs
    --------
        structure (Structure): pymatgen Structure object
                     
                            
        supercell (int):    The size supercell to use for building the 
                            bonded atomic network.
                            
        prec (float): precision of fit, in cartesian space (Angstrom)
                            
        dimension (int): number of periodic dimensions in system

        seed_index (int): starting atom for indexing                                   

    returns
    --------
        in_network (list): entries of candidate lattice vectors
    
    '''  


    # Get position of reference atom, and periodic images
    ref_spec = structure.sites[seed_index].specie
    ref_site = np.array(structure.sites[seed_index].coords)
    ref_lat  = np.array(structure.lattice.as_dict()['matrix'])/supercell

    ref_coords = [np.array(x.coords)-ref_site for x in structure.sites 
                  if x.specie ==ref_spec
                  and not magni(np.array(x.coords)-ref_site)<=prec]


    # Get the list of vectors between 
    # isolated atomic network
    
    min_area = min([magni(np.cross(x[0],x[1])) for x in 
                    list(itertools.product(ref_lat,repeat=2))
                    if magni(np.cross(x[0],x[1]))>0])
    max_area = (2**.5)*max([magni(np.cross(x[0],x[1])) for x in 
                    list(itertools.product(ref_lat,repeat=2))])

    fin_vec_list = []
    
    # Create a list of lists. Each sublist contains vectors which are all
    # parallel to each other. Each sublist is not parallel to anhy other
    # sublist
    for vec in ref_coords:
        # Check if any seen vectors are parallel
        if not np.any([(abs(np.pi-getAngle(vec,x[0][0]))<prec or 
                        abs(getAngle(vec,x[0][0]))<prec) for 
                        x in fin_vec_list]):
            fin_vec_list.append([[vec,magni(vec)]])
        else:
            i=0
            for vec_set in fin_vec_list:
                if abs(np.pi-getAngle(vec,vec_set[0][0]))<prec or \
                   abs(getAngle(vec,vec_set[0][0]))<prec:
                    fin_vec_list[i].append([vec,magni(vec)])
                    break
                i+=1
                
                
    # Get possible lattice vectors. If three vectors lie in a plane, then
    # they are grouped together as possible lattice vectors
                
    valid_vecs = []
    for i in range(len(ref_coords)):
        done=False
        for j in range(i,len(ref_coords)):
            if getAngle(ref_coords[i],ref_coords[j])>1E-2 and \
             abs(np.pi-getAngle(ref_coords[i],ref_coords[j]))>1E-2:
                cross = np.cross(ref_coords[i],ref_coords[j])
                for k in range(j,len(ref_coords)):
  
                        angs = [getAngle(ref_coords[k],cross)]
 
                        if np.all([np.pi/2-x for x in angs]):
                            done=True
                            break

            if done:
                valid_vecs.append([ref_coords[x] for x in [i,j,k]])
                
                break
 
    in_network=[]
    
    # Screen any set of vecots which exceed area criteria and reorganize
    # the data
    for v_set in valid_vecs:
        for pair in list(itertools.product(v_set,repeat=2)):
            cross_mag=magni(np.cross(pair[0],pair[1]))
            if cross_mag>= min_area and cross_mag<=max_area:
                in_network.append([pair[0],pair[1],
                                   getAngle(pair[0],pair[1])-np.pi/2,
                                   cross_mag])


    
    

    in_network = sorted(in_network, key=itemgetter(3))
    return(in_network) 
  



def getVectors(structure, supercell, prec,dim,seed_index=0):
    

    '''
    
    Helper function to get the in-plane lattice vectors.
    
        inputs
    --------
        structure (Structure): pymatgen Structure object
                     
                            
        supercell (int):    The size supercell to use for building the 
                            bonded atomic network.
                            
        prec (float): precision of fit, in cartesian space (Angstrom)
                            
        dimension (int): number of periodic dimensions in system

        seed_index (int): starting atom for indexing                                   

    returns
    --------
        in_network (list): entries of candidate lattice vectors
    
    '''  


    # Get position of reference atom, and periodic images
    ref_spec = structure.sites[seed_index].specie
    ref_site = np.array(structure.sites[seed_index].coords)
    ref_lat  = np.array(structure.lattice.as_dict()['matrix'])/supercell

    ref_coords = [[np.array(ref_site+np.dot(ref_lat.T,x)),x] for x in 
                 list(itertools.product([1,-1,0],repeat=3))]
    cluster_sites = [i for i in structure.sites]


    # Get the list of reference atoms which are in the
    # isolated atomic network
    
    in_network = []
    in_network_indexes=[]
    i=0
    test_vecs = []

    for atom in ref_coords:
        if np.any([magni(atom[0]-x.coords)<prec for x in 
                   cluster_sites if x.specie==ref_spec]):
            in_network.append(atom[0])
            in_network_indexes.append(i)
            test_vecs.append(atom[1])
        i+=1
        
    if dim==2:
        all_angs = []
    
        if len(test_vecs)>=2:
            tv = np.cross(test_vecs[0],test_vecs[1])
            for t in [x for x in test_vecs if magni(x)!=0]:
                temp_ang = getAngle(tv,t)
                if temp_ang not in all_angs:
                    all_angs.append(temp_ang)
    
        fresh_angs=[]
    
        if len(all_angs)>1:
            fresh_angs= []
            for ang in all_angs:
                if not np.any([(ang-x)<prec for x in fresh_angs]):
                    fresh_angs.append(ang)
        elif len(all_angs)>0:
            fresh_angs=[all_angs[0]]
    
        if len(fresh_angs)==1:
            new = [[np.dot(ref_lat.T,x[0]),np.dot(ref_lat.T,x[1]),
                   abs(getAngle(x[0],x[1])-np.pi/2),
                   magni(np.cross(np.dot(ref_lat.T,x[0]),
                                  np.dot(ref_lat.T,x[1])))] 
                   for x in list(itertools.combinations(
                           [v for v in test_vecs if magni(v)>0],dim))]
            return(new)
        else:
            return([])
    else:
        return([[np.dot(ref_lat.T,v),np.dot(ref_lat.T,v),np.pi/4,
                 magni(v),np.dot(ref_lat.T,v)] for v in test_vecs 
                if magni(v)>0])
    return(in_network)

def getNewLattice(entry,dim,prec=1E-4,seed_index=0,supercell=4,c_mag=60,
                  attempt_rough=True):
    
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
        
        attempt_rough (bool): If no in-plane lattice vectors are found which
                              allow a single unit cell representation for 
                              the monolayer, search for any combination of 
                              lattice vectors which fit the crystal structure.
                              Typically only necessary for cleaved structures 
                              where the surface cleaved needs more than one
                              unit cell of the precursor to be formed 
                              (EX: in-plane vectors needed are [0,1,0] and 
                               [2,.5,1], dotted with the overall lattice 
                               vectors)

        
    returns
    --------
        latt_attempt (list, 3X3): New lattice for "structure"
        
    """   
    
    
    
    # Get network of atoms within the supercell

    # Reset structure object, as "getStructureType" turned it into a supercell
    

    structure = entry[0]
    
    # Shift atoms so that the refernce atom is at the center of the cell
    structure.translate_sites(indices=range(structure.num_sites),
                            vector=-1*structure.sites[seed_index].frac_coords+
                            [0.99,0.99,0.99],frac_coords=True)
                             
    structure = reduceCoords(structure)
    ogStructure = copy.deepcopy(structure)
    
    # Replace structure used for TSA with centered structure
    entry[0]=structure
    batch = getStructureType(entry,seed_index=seed_index,
                             supercell=supercell,return_SS=True)

    structure.make_supercell(supercell)
    cluster_sites = [structure.sites[n] for n in batch[4]]
    
    
    

    # Get structure object containing a single atomic network
    
    
    ogStructure = Structure.from_sites(cluster_sites)
    ogStructure = reduceCoords(ogStructure)
    in_network = getVectors(ogStructure,supercell,prec,2)
    fit_fracs=[]

    # Get the new lattice matrix
    latt_attempt,fit_fracs= genLattice(ogStructure,
                                       in_network,dim,supercell,c_mag=c_mag)
    
    
    # If a valid lattice matrix isn't found, it is likely due to an issue
    # with original lattice vectors not being orthogonal. 
    # Shift atomic positions a bit until this problem doesn't occur
    o_shift=0.05
    o_shift_count=float(o_shift)
    
    # While no matching set of lattice vectors is found,
    # slightly shift atomic positions in an 
    # attempt to capture an accurate fit. 
    while len(latt_attempt)==0 and o_shift_count<1:

        ogStructure.translate_sites(indices=range(ogStructure.num_sites),
                                    vector=[o_shift,o_shift,o_shift],
                                    frac_coords=True)

        o_shift_count+=o_shift
        in_network = getVectors(ogStructure,supercell,prec,dim)

        latt_attempt,fit_fracs = genLattice(ogStructure,
                                    in_network,dim,supercell,c_mag=c_mag)


    # If no match was found and attempt_rough is True, attempt again with
    # a more complete list of possible lattice vectors
    if len(latt_attempt)==0 and attempt_rough:
        o_shift=0.05
        o_shift_count=float(o_shift)
        while len(latt_attempt)==0 and o_shift_count<1:
            ogStructure.translate_sites(indices=range(ogStructure.num_sites),
                                        vector=[o_shift,o_shift,o_shift],
                                        frac_coords=True)
            
            
            o_shift_count+=o_shift
            in_network = getVectorsRough(ogStructure,supercell,prec,dim)
            latt_attempt,fit_fracs = genLattice(ogStructure,
                                        in_network,dim,supercell,c_mag=c_mag)

    
    return(latt_attempt,fit_fracs)



def genLattice(structure,in_network,dim,supercell,prec=1E-4,
               seed_index=0,c_mag=60,y_dist=-1):
    
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
    
    if y_dist==-1:
        y_dist=c_mag/3
    
    new = [x for x in in_network if abs(x[2])<np.pi/2]
    return_structure=False
    mat = np.array(structure.lattice.as_dict()['matrix'])
    coords = np.array([np.dot(mat.T,x.frac_coords%1) for x in structure.sites])
    specs = structure.species
    ref_ele_d = getUniqueCount(specs)
    for i in ref_ele_d:
        ref_ele_d[i]/=(supercell**dim)
    coords = coords-coords[seed_index]
    




    for lat_vectors in sorted(new,key=itemgetter(3)):

        # Create lattice matrix to fit atomic coordinates against
        # In 2D
        if dim==2:
            new_c = np.cross(lat_vectors[0],lat_vectors[1])
            scale_c = c_mag/magni(new_c)

            latt_attempt = np.array([lat_vectors[0],lat_vectors[1],\
                                     new_c*scale_c])
                
        # In 1D
        elif dim==1:
            unitV = lat_vectors[0]/magni(lat_vectors[0])
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

        if len([x for x in np.array(new_fracs).T if 
                np.all([(y>=0 and y<1) for y in np.around(x[:dim],3)]) and
                np.all([(y>=-y_dist/c_mag and y<y_dist/c_mag) for 
                    y in np.around(x[dim:],3)])])==len(new_fracs[0])/supercell**dim:
       
                fit_fracs=[]
                new_fracs_t = np.around(new_fracs.T,6)
                for i in range(len(new_fracs[0])):
                    if np.all([(y>=0 and y<1) for y in np.around(new_fracs_t[i][:dim],3)]) \
                                  and np.all([(y>=-y_dist/c_mag and y<y_dist/c_mag) 
                                              for y in np.around(new_fracs_t[i][dim:],3)]):
                                      fit_fracs.append([new_fracs_t[i],specs[i]])
                fit_fracs = np.array(fit_fracs).T
                new_ele_d = getUniqueCount(fit_fracs[1])
                unequal=False
                for k in new_ele_d:
                    if new_ele_d[k]!=ref_ele_d[k]:
                        unequal=True

                        break
                if not unequal:

                    return_structure=True
                    break



    # If match found
    if return_structure:
        return(np.array(latt_attempt),fit_fracs)
    # If no match found
    else:
        return([],[])



def alignMono(entry,prec=1E-4,seed_index=0,supercell=2,
              c_mag=50,dist_from_plane=3):
    
    """
    Align a 2D material such that the 'c' vector is perpendicular to the
    in-plane lattice vectors. Draws lattice vectors between known atoms.
    
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
        list1 (list): -species associated with each site
                     -fractional coordinates in new lattice
                     -new lattice (a,b,c)
    
        list2 (list): -species associated with each site
                     -fractional coordinates in new lattice
                     -new lattice (a,b,c)
        
    """




    new_latt,fit_fracs_both= getNewLattice(entry,dim=2,prec=prec,
                                           seed_index=seed_index,
                                           supercell=supercell,c_mag=c_mag)

    fit_fracs = np.array([np.array(x)+[0,0,.5] for x in fit_fracs_both[0]])
    final_sites = np.dot(new_latt.T,fit_fracs.T).T
    # Create new lattice matricies
    lat1 = np.array([new_latt[0],new_latt[1],new_latt[2]])
    lat2 = np.array([new_latt[1],new_latt[0],new_latt[2]])

    # Generate atomic fractions
    new_fracs1 = np.linalg.solve(lat1.T,np.array(final_sites).T).T
    new_fracs2 = np.linalg.solve(lat2.T,np.array(final_sites).T).T
    species = fit_fracs_both[1]
    return([species,new_fracs1,lat1],[species,new_fracs2,lat2])



def alignMonoPlane(entry,prec=1E-4,seed_index=0,supercell=2,
              c_mag=50,dist_from_plane=3):
    
    """
    Align a 2D material such that the 'c' vector is perpendicular to the
    in-plane lattice vectors. Outdated in accuracy and applicability by
    the "alignMono" function.
    
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


    new_latt,fit_fracs= getNewLattice(entry,dim=2,prec=prec,seed_index=seed_index,
                          supercell=supercell,c_mag=c_mag)


    

    # Identify plane to translate atoms towards

    plane = Plane(Point3D(s.sites[seed_index].coords),
                  normal_vector=new_latt[2])
    
    # Create list of translationss
    trans = list(itertools.product([1,-1,0],repeat=3))

    lat = np.array(s.lattice.as_dict()['matrix'])
    final_sites = []
    i=0
   
    # Ensure that the atoms are nearby each other
    for site in [x.coords for x in s.sites]:
        point = Point3D(site)
        if 1==1:

            news = []
            
            # translate atomic sites to see which position is closest to plane
            for t in trans:
                point = Point3D(site+np.dot(np.transpose(lat),t))
                news.append([float(plane.distance(point)),t])
            news.sort(key = lambda x:x[0])
            for new in news:
                if not np.any([magni((site+np.dot(np.transpose(lat),new[1]))-x)<=prec for x in final_sites]):
                    final_sites.append(site+
                                       np.dot(np.transpose(lat),new[1]))
                    break
        i+=1
        
    # Create new lattice matricies
    lat1 = np.array([new_latt[0],new_latt[1],new_latt[2]])
    lat2 = np.array([new_latt[1],new_latt[0],new_latt[2]])

    # Generate atomic fractions
    new_fracs1 = np.linalg.solve(lat1.T,np.array(final_sites).T).T
    new_fracs2 = np.linalg.solve(lat2.T,np.array(final_sites).T).T

    species=fit_fracs[1]

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

    v1,v2,perp=new_latt
    new_latt=np.array(new_latt)
    line1 = Line3D(new_latt[0],[0,0,0])
    line = line1.parallel_line(new_struct.sites[0].coords)
    trans = list(itertools.product([1,-1,0],repeat=3))

    lat = np.array(new_struct.lattice.as_dict()['matrix'])
    final_sites = []
    i=0
    i=-1
    for site in [x.coords for x in new_struct.sites]:

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
        ang = abs(getAngle(new_latt[0],new_latt[1]))
        if a>b:
            frac_coords = [[x[1],x[0],x[2]] for x in frac_coords]

            a,b=b,a
        if abs((ang-np.pi/2)/ang)<1E-4:
            new_latt = [[a,0,0],
                        [0,b,0],
                        [0,0,c]]
        else:
            new_latt=[[a*np.cos(ang/2),a*np.sin(ang/2),0],
                      [b*np.cos(ang/2),-1*b*np.sin(ang/2),0],
                      [0,0,c]]
    elif dim==1:    
        new_latt = [[a,0,0],[0,b,0],[0,0,c]]
    i=0
    
    new_sites = []
    for site in frac_coords:
        new_coords = np.dot(np.array(new_latt).T,np.array(site).T).T
        p = PeriodicSite(species=Element(specs[i]),
                         lattice = Lattice(new_latt),
                         coords=new_coords,
                         coords_are_cartesian=True)
        new_sites.append(p)
        i+=1

    new_struct = Structure.from_sites(new_sites)

    return(new_struct)



def reduceScale(structure,scale,dim,prec=1E-2):
    
    """
    Attempt to reduce the structure from a supercell
    
    inputs
    --------
        structure (Structure):  Pymatgen structure object to reduce
        
        scale (int or list (3x1)): How much to scale the periodic lattice vectors by.
                               If int, scales all periodic lattice vectors by
                               that int. If list, scales lattice vectors by the
                               entries of the list (x,y,z)
        
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
    if type(scale)==type(1) or type(scale)==type(1.0):
        lat[0]*=scale
        if dim>=2:
            lat[1]*=scale
        if dim==3:
            lat[2]*=scale
    else:
        lat[0]*=scale[0]
        lat[1]*=scale[1]
        lat[2]*=scale[2]        

    fracs = np.around(abs(np.array(np.linalg.solve(lat.T,np.array(cart).T).T)),4)

    specs   = []
    u_cart   = []
    i=0
    n_fracs = []
    for frac in fracs:
        if not np.any([magni(np.around(frac,4)%1-np.around(x,4)%1)<prec for x in n_fracs]):
                n_fracs.append(frac%1)
        #if frac[0]<1 and frac[1]<1 and frac[2]<1:
            
                specs.append(structure.species[i])
                u_cart.append(cart[i])
        i+=1
    
    new_sites = []
    i=0
    for site in u_cart:
        p = PeriodicSite(lattice = Lattice(lat),coords=site,
                         coords_are_cartesian=True,
                         species=specs[i])
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

