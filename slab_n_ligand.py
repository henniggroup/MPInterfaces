import sys
import numpy as np
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.operations import SymmOp

def create_slab_n_adsorb(strt, hkl, min_thick, min_vac, supercell, mol, shift):
    slb = SlabGenerator(strt, hkl, min_thick, min_vac).get_slab()   
    slb.to(fmt='poscar', filename='POSCAR_primitive_slab.vasp')
    #print slb.surface_area
    #create supercell of the slab
    slb.make_supercell(supercell)
    slb.to(fmt='poscar', filename='POSCAR_slab_supercell.vasp')
    #get all the top atom coords
    cart_coords = slb.cart_coords
    max_dist = np.max(slb.distance_matrix[0,:])
    #print 'max_dist', max_dist
    #get the list of top atoms
    top_atoms = []
    for j in range(len(cart_coords[:,0])):
        if j not in top_atoms:
            [top_atoms.append(i) for i in range(len(cart_coords[:,0])) if (max_dist-slb.distance_matrix[j,i]) < 1e-6]
    #top_atom_coords = [cart_coords[top_atom] for top_atom in top_atoms]
    #print 'top_atoms, top_atom_coords', top_atoms, top_atom_coords
    #set adatoms_coords wrt the top atoms in Angstrom
    adatoms_coords = []
    for top_atom in top_atoms:
        adatoms_coords.append(mol.cart_coords + cart_coords[top_atom] + shift)

    adatoms_coords = np.array(adatoms_coords) #3d numpy array
    #print adatoms_coords
    #extend the slab structure with the adsorbant atoms
    for j in range(len(top_atoms)):
        [slb.append(adatoms[i], adatoms_coords[j,i,:], coords_are_cartesian=True) for i in range(len(adatoms))]
    #print slb.cart_coords
    
    return slb
    

#test
if __name__=='__main__':
    #latt const from materialsproject
    a0 = 3.62
    #conventional cell and h,k,l wrt that
    latt = Lattice.cubic(a0)
    species = ["Cu", "Cu", "Cu", "Cu"]
    positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    
    hkl = [1,1,1]
    min_thick = 9 #Ang
    min_vac = 12 #Ang
    
    strt = Structure(latt, species, positions)
    supercell = [2,2,1]
    
    # PBEPBE/aug-cc-pVQZ , from cccbdb
    #0 0 0.119079
    #0 0.7648 -0.476318
    #0 -0.7648 -0.476318
    #adsorb species 'O' on atom of index 3 shifted by 1.0 A
    #slb.add_adsorbate_atom([3], 'O', 1.0)
    adatoms = ['O','H', 'H']
    adatoms_coords = np.array([ [0,0,0], [0, 0.77, 0.60], [0, -0.77, 0.60]])
    mol = Molecule(adatoms, adatoms_coords)
    #print mol
    #rotate along the specified axis
    symop = SymmOp.from_axis_angle_and_translation([1,0,0], -45)
    mol.apply_operation(symop)
    shift = [0, 0, 2]
    #print mol

    #create a slab supercell of miller index hkl, minimum thickness min_thick and a
    #minimum vacuum space of min_vac
    #adsorb a molecule(shifted by shift from the slab top surface) on
    #NOTE: the molecule is adsorbed on all atoms of the top surface of the slab
    #TODO: add surface coverage(ligands per nm^2)
    slb = create_slab_n_adsorb(strt, hkl, min_thick, min_vac, supercell, mol, shift)
    
    #sort and write to poscar
    slb.sort()
    slb.to(fmt='poscar', filename='POSCAR_final.vasp')
    

