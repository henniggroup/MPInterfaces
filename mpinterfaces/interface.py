import sys
import numpy as np
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.core.operations import SymmOp

class Interface(Slab):
    """
    Interface = slab + ligand + environment(solvent)
    """
    def __init__(self, strt, hkl=[1,1,1], min_thick=10, min_vac=10, supercell=[1,1,1],
                 ligand=None, displacement=1.0, sface_cvrg=None, solvent=None, start_from_slab=False,
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
        self.sface_cvrg = sface_cvrg
        
        Slab.__init__(self, self.strt.lattice, self.strt.species_and_occu, self.strt.frac_coords, miller_index=self.hkl, oriented_unit_cell=None, shift=None, scale_factor=None, validate_proximity=validate_proximity,
                           to_unit_cell=to_unit_cell, coords_are_cartesian=coords_are_cartesian,
                           site_properties=self.strt.site_properties, energy=None )

    def set_top_atoms(self):
        """
        get all the top atom coords
        """
        cart_coords = self.strt.cart_coords
        n_atoms = len(cart_coords[:,0])
        max_dist = np.max(self.strt.distance_matrix.reshape(n_atoms*n_atoms, 1))
        self.top_atoms = []
        for j in range(n_atoms):
            if j not in self.top_atoms:
                [self.top_atoms.append(i) for i in range(n_atoms) if np.abs(max_dist - self.strt.distance_matrix[j,i]) < 1e-6]
        self.top_atoms = np.unique(self.top_atoms)

    def enforce_surface_cvrg(self):
        m = self.strt.lattice.matrix
        surface_area = np.linalg.norm(np.cross(m[0], m[1]))
        #n_top_atoms =  len(self.top_atoms)
        #print 'number of top atoms ', n_top_atoms
        if self.sface_cvrg:
           for i in range(1, 10):
                for j in range(1, 10):
                    #for nlig in range(1, i*j*n_top_atoms):
                        
                        surface_area = np.linalg.norm(np.cross(i*m[0], j*m[1]))
                        if (int(surface_area * self.sface_cvrg) == 1):
                                print 'supercell = ', i, j
                                print 'feasible covergae = ', 1./surface_area, ' requested = ', self.sface_cvrg
                                self.strt.make_supercell([i,j,1])
                                return


    def adsorb_on(self, site_indices):
        mov_vec = self.strt.normal * self.displacement
#        for i in range(len(self.strt.cart_coords)):
#            if np.linalg.norm(self.strt.cart_coords[i] - mov_vec) <= 3.0:
#                mov_vec = -mov_vec
#                break
        #print self.strt.cart_coords[site_indices[0]]
        #sys.exit()
        
        adsorbed_ligands_coords = np.array([self.ligand.cart_coords - self.strt.cart_coords[sindex] + mov_vec for sindex in site_indices ]) #3d numpy array
        #extend the slab structure with the adsorbant atoms
        num_atoms = len(self.ligand.species_and_occu)
        for j in range(len(site_indices)):
            [self.strt.append(self.ligand.species_and_occu[i], adsorbed_ligands_coords[j,i,:], coords_are_cartesian=True) for i in range(num_atoms)]
        
                
    def create_interface(self):
        if not self.start_from_slab:
            self.strt = SlabGenerator(self.strt, self.hkl, self.min_thick, self.min_vac).get_slab()
        self.strt.to(fmt='poscar', filename='POSCAR_primitive_slab.vasp')
        #self.strt.to(fmt='poscar', filename='POSCAR_slab_supercell.vasp')
        #self.set_top_atoms()
        #self.adsorb_sites = self.top_atoms
        print self.strt
        self.enforce_surface_cvrg()
        print self.strt
        self.set_top_atoms()
        self.adsorb_sites = self.top_atoms
        print self.adsorb_sites
        sys.exit()
        #adsorb one ligand on one of the top atoms
        self.adsorb_on([self.adsorb_sites[0]])        
        #sys.exit()
        #print 'surface coverage ', 1./self.strt.surface_area


        #create supercell of the slab to enforce the given surface coverage
        self.enforce_surface_cvrg()


    @staticmethod
    def from_file(fname):
        pass

    def write_to_file(self,fmt,fname):
        self.strt.sort()
        self.strt.to(fmt=fmt, filename=fname)


        
class Ligand(Molecule):
    """
    Construct ligand from  molecules
    
    """
    def __init__(self, mols, cm_dist=[], angle={}, link={}, remove=[],
                 charge=0, spin_multiplicity=None,
                 validate_proximity=False):
        Molecule.__init__(self, mols[0].species_and_occu, mols[0].cart_coords,
                          charge=charge, spin_multiplicity=spin_multiplicity,
                          validate_proximity=validate_proximity,
                          site_properties=mols[0].site_properties)
        self._sites = list(self._sites)
        self.mols = mols
        self.cm_dist = cm_dist
        self.angle = angle
        self.link = link
        self.remove = remove
        

    def get_perp_vec(self, vec1, vec2):
        #if the vectors are parllel, then perp_vec = (0, -z, y)
        if np.abs(np.dot(vec1, vec2) - np.linalg.norm(vec1)**2 ) < 1e-6:
            perp_vec = np.array([0, -vec1[2], vec1[1]])
        else:
            perp_vec = np.cross(vec1, vec2)
        return perp_vec                        


    
    def get_distance_matrix(self, mol):
        nsites = len(mol.sites)
        return np.array([mol.get_distance(i,j) for i in range(nsites) for j in range(nsites)]).reshape(nsites, nsites)



    def create_ligand(self):
        #get the start and end indices to define the vector that defines the molecule
        self.vec_indices = []
        for mol in self.mols:
            nsites = len(mol.sites)
            d_mat = self.get_distance_matrix(mol)
            max_dist = np.max(d_mat.reshape(nsites*nsites, 1))
            temp = []
            for i in range(nsites):
                if i not in temp:
                    [temp.append([i,j]) for j in range(nsites) if np.abs(max_dist-d_mat[i,j]) < 1e-6]
            self.vec_indices.append(temp[0])

        print self.vec_indices
        #vectors pointing from start index atom to the farthest atom
        self.mol_vecs = []
        for mol,vind in enumerate(self.vec_indices):
            self.mol_vecs.append(self.mols[mol].cart_coords[vind[1]] - self.mols[mol].cart_coords[vind[0]])
    
        #move center of masses
        new_mol = self.mols[0]
        #first movemen tin the x direction
        mov_vec = np.array([1,0,0])
        for i in range(len(self.mols)-1):
            cm1 = new_mol.center_of_mass
            new_cm = new_mol.center_of_mass        
            cm2 = self.mols[i+1].center_of_mass
            new_cm = new_cm + self.cm_dist[i] * mov_vec #+ np.random.rand(1,3)
            mov_vec = self.get_perp_vec(self.mol_vecs[i], mov_vec)
            mov_vec = mov_vec / np.linalg.norm(mov_vec)
            new_coords = self.mols[i+1].cart_coords + new_cm
            self.mols[i+1] = Molecule(self.mols[i+1].species_and_occu, new_coords,
                          charge=self.mols[i+1]._charge,
                          spin_multiplicity=self.mols[i+1]._spin_multiplicity,
                          site_properties=self.mols[i+1].site_properties)
            new_mol = Molecule.from_sites(self.mols[i].sites + self.mols[i+1].sites, validate_proximity=True) 

            #rotate the molecules around an axis that is
            #perpendicular to the molecular axes
            if self.angle:
                for mol in range(len(self.mols)):
                    for ind_key, rot in self.angle[str(mol)].items():
                        #print 'mol, ind_key, rot ', mol, ind_key, rot
                        perp_vec = np.cross(self.mol_vecs[int(ind_key)], self.mol_vecs[mol])
                        #if the vectors are parllel, then perp_vec = (-y, x, 0)
                        if np.abs(np.dot(self.mol_vecs[int(ind_key)], self.mol_vecs[mol]) - np.linalg.norm(self.mol_vecs[mol])**2 ) < 1e-6:
                            perp_vec = np.array([-self.mol_vecs[mol][1], self.mol_vecs[mol][0], 0])
                        org_pt = self.vec_indices[mol][0]
                        op = SymmOp.from_origin_axis_angle(self.mols[mol].cart_coords[org_pt], axis=perp_vec, angle=rot)
                        self.mols[mol].apply_operation(op)
            #connect the molecules together
            #in this example: connect mol2 to mol0 and mol1, the
            #coordinates of mol2 is changed 
            new_coords = np.array([0,0,0]) 
            displacement = np.array([0,0,0])
            if self.link:
                for mol in range(len(self.mols)):
                    new_coords = copy.deepcopy(self.mols[mol].cart_coords)
                    if link[str(mol)]:
                        for ind_key, conn in self.link[str(mol)].items():
                            ind = int(ind_key)
                            print 'connection list for atom of index ',ind ,' of molecule ', mol, ' : ', conn
                            coord = np.array([0,0,0])
                            #if connecting the molecule mol to only one atom of
                            #just one another molecule
                            #then move the atom close to the atom in mol and
                            #shift the whole molecule too
                            non_neg = np.extract(np.array(conn)>0, conn)
                            if len(non_neg) == 1 and len(link[str(mol)]) == 1:
                                for j,k in enumerate(conn):
                                    coord = self.mols[j].cart_coords[non_neg[0]] + np.random.rand(1,3) + 1.0
                                displacement = coord - self.mols[mol].cart_coords[ind]
                                #connect the specified atoms of mol to the atoms of
                                #other molecules in the list
                                #connection means putting the atomm  of the mol at
                                #a position that is the average of the position of
                                #the atoms of the molecules given in the list
                            else:
                                for j,k in enumerate(conn):
                                    if k>=0:
                                        coord = coord + self.mols[j].cart_coords[k] 
                                coord = coord / len(conn)
                            new_coords[ind, 0] = coord[0]
                            new_coords[ind, 1] = coord[1]
                            new_coords[ind, 2] = coord[2]
                        new_coords = new_coords + displacement
                        self.mols[mol] =  Molecule(self.mols[mol].species_and_occu,
                                   new_coords, charge=self.mols[mol]._charge,
                                   spin_multiplicity=self.mols[mol]._spin_multiplicity,
                                   site_properties=self.mols[mol].site_properties)
    
        #remove atoms from the molecules
        for i in range(len(self.mols)):
            if self.remove[i]:
                self.mols[i].remove_sites(self.remove[i])

        #combine the sites of all the molecules
        combine_mol_sites = self.mols[0].sites
        for j in range(1, len(self.mols)):
            combine_mol_sites = combine_mol_sites + self.mols[j].sites
        
        #create the ligand
        self._sites = combine_mol_sites


#test
if __name__=='__main__':
    ##########################
    #lead acetate ligand
    ##########################
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
    lead_acetate = Ligand(mols, cm_dist, angle=angle, link={}, remove=remove)
    lead_acetate.create_ligand()
    lead_acetate.to('xyz', 'lead_acetate.xyz')

    #put the ligand in a box
    boxed_lead_acetate = lead_acetate.get_boxed_structure(13, 13, 13)  
    boxed_lead_acetate.to(fmt= "poscar", filename= "POSCAR_diacetate_boxed.vasp")

    ######################################
    # create a slab corresponding to an
    #  interface with h2o lignads
    # attached to Cu 111 surface    
    # Cu bulk latt const from materialsproject
    #####################################
    a0 = 3.62
    latt = Lattice.cubic(a0)
    species = ["Cu", "Cu", "Cu", "Cu"]
    positions = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]

    ########################################
    # H2O ligand
    ########################################        
    # PBEPBE/aug-cc-pVQZ , from cccbdb
    #0 0 0.119079
    #0 0.7648 -0.476318
    #0 -0.7648 -0.476318
    #adsorb species 'O' on atom of index 3 shifted by 1.0 A
    #slb.add_adsorbate_atom([3], 'O', 1.0)
    adatoms = ['O','H', 'H']
    adatoms_coords = np.array([ [0,0,0], [0, 0.77, 0.60], [0, -0.77, 0.60]])
    mol = Molecule(adatoms, adatoms_coords)
    symop = SymmOp.from_axis_angle_and_translation([1,0,0], -45)
    mol.apply_operation(symop)
    h2o = Ligand([mol])

    ########################################
    # Interface = Slab + ligand + solvent
    #######################################
    hkl = [1,1,1]
    min_thick = 9 #Ang
    min_vac = 12 #Ang    
    strt = Structure(latt, species, positions)
    supercell = [2,2,1]
    #ligand displacement from the slab surface along the surface normal
    displacement = 2.0 

    #create a slab supercell of miller index hkl, minimum thickness min_thick and a
    #minimum vacuum space of min_vac
    #adsorb a molecule(shifted by shift from the slab top surface) on
    #NOTE: the molecule is adsorbed on all atoms of the top surface of the slab
    #TODO: add surface coverage(ligands per nm^2)
    iface = Interface(strt, hkl=hkl, min_thick=min_thick, min_vac=min_vac, supercell=supercell, ligand=h2o, displacement=displacement, sface_cvrg=0.1) # 1 lig/nm^2 = 0.01 lig/ang^2
    iface.create_interface()
    iface.write_to_file('poscar', 'POSCAR_interface.vasp')



    

