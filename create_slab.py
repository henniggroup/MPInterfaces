import sys
import numpy as np
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.operations import SymmOp


if __name__=='__main__':
    #latt const from materialsproject
    a0 = 3.62
    #conventional cell and h,k,l wrt that
    latt = Lattice.cubic(a0)
    species = ["Cu", "Cu", "Cu", "Cu"]
    positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    
    hkl = [1,1,1] # miller indices
    min_thick = 9 # in Angstrom
    min_vac = 10 # in Angsstrom
    
    strt = Structure(latt, species, positions)
    slb = SlabGenerator(strt, hkl, min_thick, min_vac).get_slab()   
    slb.to(fmt='poscar', filename='POSCAR_cu_111.vasp')
    print slb.cart_coords

    # PBEPBE/aug-cc-pVQZ , from cccbdb database
    # 0 0 0.119079
    # 0 0.7648 -0.476318
    # 0 -0.7648 -0.476318
    #adsorb species 'O' on atom of index 3 shifted by 1.0 A
    #slb.add_adsorbate_atom([3], 'O', 1.0)
    
    adatoms = ['O','H', 'H']
    adatoms_coords = np.array([ [0,0,0], [0, 0.77, 0.60], [0, -0.77, 0.60]])
    mol = Molecule(adatoms, adatoms_coords)
    print mol

    #rotate along the specified axis
    symop = SymmOp.from_axis_angle_and_translation([1,0,0], -45)
    mol.apply_operation(symop)
    print mol

    #get the top atom coords
    cart_coords = slb.cart_coords
    top_atom = np.argmax(slb.distance_matrix[0,:])
    top_atom_coords = cart_coords[top_atom]

    #set adatoms_coords wrt the top atom in Angstom
    shift = [0, 0, 2] #shift vecor
    adatoms_coords = mol.cart_coords + top_atom_coords + shift

    #add the molecule to the slab structure
    [slb.append(adatoms[i], adatoms_coords[i,:], coords_are_cartesian=True) for i in range(len(adatoms))]
    print slb.cart_coords
    
    slb.sort()
    slb.to(fmt='poscar', filename='POSCAR_cu_111_ad.vasp')
