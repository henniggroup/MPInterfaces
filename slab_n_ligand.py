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
    def __init__(self, mols, cm_dist, angle={}, link={}, remove=[]):
        self.mols = mols
        Structure.__init__(self, self.mols[0].lattice, self.mols[0].species_and_occu, self.mols[0].frac_coords,
                      validate_proximity=validate_proximity,
                      to_unit_cell=to_unit_cell, coords_are_cartesian=coords_are_cartesian,
                      site_properties=self.mols[0].site_properties)

    def get_perp_vec(vec1, vec2):
            #if the vectors are parllel, then perp_vec = (0, -z, y)
        if np.abs(np.dot(vec1, vec2) - np.linalg.norm(vec1)**2 ) < 1e-6:
            perp_vec = np.array([0, -vec1[2], vec1[1]])
        else:
            perp_vec = np.cross(vec1, vec2)
        return perp_vec                        

    def get_distance_matrix(mol):
        nsites = len(mol.sites)
        return np.array([mol.get_distance(i,j) for i in range(nsites) for j in range(nsites)]).reshape(nsites, nsites)



    def create_ligand(self):
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

    

