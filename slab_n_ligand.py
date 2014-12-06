import sys
import numpy as np
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.core.operations import SymmOp

class Interface(Slab):
    def __init__(self, strt, hkl=[1,1,1], min_thick=10, min_vac=10, supercell=[1,1,1],
                 ligand=None, displacement=1.0, solvent=None, start_from_slab=False,
                 validate_proximity=False, to_unit_cell=False, coords_are_cartesian=False):
        self.strt = strt
        self.hkl = hkl
        self.min_thick = min_thick
        self.min_vac = min_vac
        self.supercell = supercell
        self.ligand = ligand
        self.displacement = displacement
        self.solvent = solvent
        self.start_from_slab = start_from_slab
        
        Slab.__init__(self, self.strt.lattice, self.strt.species_and_occu, self.strt.frac_coords, miller_index=self.hkl, oriented_unit_cell=None, shift=None, scale_factor=None, validate_proximity=validate_proximity,
                           to_unit_cell=to_unit_cell, coords_are_cartesian=coords_are_cartesian,
                           site_properties=self.strt.site_properties, energy=None )

    def set_top_atoms(self):
        """
        get all the top atom coords
        """
        cart_coords = self.strt.cart_coords
        max_dist = np.max(self.strt.distance_matrix[0,:])
        self.top_atoms = []
        for j in range(len(cart_coords[:,0])):
            if j not in self.top_atoms:
                [self.top_atoms.append(i) for i in range(len(cart_coords[:,0])) if (max_dist - self.strt.distance_matrix[j,i]) < 1e-6]

                
    def create_interface(self):
        if not self.start_from_slab:
            self.strt = SlabGenerator(self.strt, self.hkl, self.min_thick, self.min_vac).get_slab()
        self.strt.to(fmt='poscar', filename='POSCAR_primitive_slab.vasp')
        #create supercell of the slab
        self.strt.make_supercell(supercell)
        self.strt.to(fmt='poscar', filename='POSCAR_slab_supercell.vasp')
        normal =  self.strt.normal
        self.set_top_atoms()
        #set adatoms_coords wrt the top atoms in Angstrom
        adsorbed_ligands_coords = np.array([self.ligand.cart_coords + self.strt.cart_coords[top_atom] + normal * self.displacement for top_atom in self.top_atoms ]) #3d numpy array
        #extend the slab structure with the adsorbant atoms
        num_atoms = len(self.ligand.species_and_occu)
        for j in range(len(self.top_atoms)):
            [self.strt.append(self.ligand.species_and_occu[i], adsorbed_ligands_coords[j,i,:], coords_are_cartesian=True) for i in range(num_atoms)]


    @staticmethod
    def from_file(fname):
        pass

    def write_to_file(self,fmt,fname):
        self.strt.sort()
        self.strt.to(fmt=fmt, filename=fname)


        
class Ligand(Structure):
    def __init__(self):
        pass

    def create_ligand(self):
        pass


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
    displacement = 2.0 #[0, 0, 2]
    #print mol

    #create a slab supercell of miller index hkl, minimum thickness min_thick and a
    #minimum vacuum space of min_vac
    #adsorb a molecule(shifted by shift from the slab top surface) on
    #NOTE: the molecule is adsorbed on all atoms of the top surface of the slab
    #TODO: add surface coverage(ligands per nm^2)

    iface = Interface(strt, hkl=hkl, min_thick=min_thick, min_vac=min_vac, supercell=supercell, ligand=mol, displacement=displacement)
    iface.create_interface()
    iface.write_to_file('poscar', 'POSCAR_final.vasp')

    

