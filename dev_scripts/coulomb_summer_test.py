"""
Mind: Specify  path for your bader executable in main with the needed options
"""

from pymatgen import Structure, Molecule
from pymatgen.analysis.energy_models import EnergyModel
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.io.vaspio import VolumetricData, Chgcar, Potcar, Locpot
import os
import numpy as np
from mpinterfaces import *
from pymatgen.io.vaspio import Poscar
from ase.io.vasp import write_vasp
from pymatgen.io.aseio import AseAtomsAdaptor
from ase.calculators.calculator import Calculator




class ColoumbEnergy(EnergyModel):
    def get_energy(self, structure):
        """
        use the atomic numbers of the structures's elements to compute 
        the coloumbic energy of the structure
        """
        energy = 0
        a = -1
        b = -1
        for i, sitei in enumerate(structure):#enumerate(structure.as_dict()['sites']):
            for j, sitej in enumerate(structure):#enumerate(structure.as_dict()['sites']):
                if i != j:
                    dij = structure.get_distance(i,j)
                    Zi=sitei.species_and_occu.items()[0][0].Z
                    Zj=sitej.species_and_occu.items()[0][0].Z
                    energy +=  0.5 * Zi*Zj/dij
        return energy
     
    def get_bader_coulomb_energy(self, structure):
	"""
	use bader charges on atoms to compute the coloumbic energy
	"""
        energy = 0
        force = np.zeros((len(structure), 3))
        a = -1
        b = -1
        for i, sitei in enumerate(structure.as_dict()['sites']):
            for j, sitej in enumerate(structure.as_dict()['sites']):
                if i != j:
                    dij = structure.get_distance(i,j)
		    d_vec = structure.frac_coords[i] - structure.frac_coords[j]
                    Zi = sitei['species'][0]['oxidation_state']
                    Zj = sitej['species'][0]['oxidation_state']
                    energy +=  0.5 * Zi*Zj/dij
	    	    force[i][0] += Zi*Zj/(dij**2) * (d_vec[0]/dij)
	    	    force[i][1] += Zi*Zj/(dij**2) * (d_vec[1]/dij)
		    force[i][2] += Zi*Zj/(dij**2) * (d_vec[2]/dij)
	    print force[i]
				 #to work on definition of forces 
        print np.sum(force[:,0]), np.sum(force[:,1]), np.sum(force[:,2]) #total force on cell in x, y, z ? 
	return energy

    def get_ewald_sum(self,structure):
        e = EwaldSummation(structure, real_space_cut=None,
                           recip_space_cut=None,
                           eta=None,
                           acc_factor=8.0)
        return e.total_energy
 
    def as_dict(self):
        return {"version": __version__,
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "init_args": {"name":"coloumb"} }



class Bader_Analysis(object):
    """
    """
    def __init__(self, acf_path= './ACF.dat', chgcar_filename= "./CHGCAR", aecar0 = "./AECCAR0", aecar2 = "./AECCAR2", potcar_filename="./POTCAR", bader_path= "./bader CHGCAR -ref Charge_sum"):
        print "Reading CHGCAR"
        self.chgcar = Chgcar.from_file(chgcar_filename)
	##uncomment if you have to run from scratch##

#	self.contcar = Structure.from_file(contcar_filename)
#	print "Reading AECCAR0"
#       Vol_obj_1 = Chgcar.from_file(aecar0)
#	print "Reading AECCAR2"
#       Vol_obj_2 = Chgcar.from_file(aecar2)
#	print "Summing"
#       Vol_obj_sum = Vol_obj_1.linear_add(Vol_obj_2)
#	print "Writing Combined Sum"
#       Vol_obj_sum.write_file("./Charge_sum")
#	self.exe = bader_path
#	os.system(self.exe)
        self.potcar = Potcar.from_file(potcar_filename)             if potcar_filename is not None else None
        self.natoms = self.chgcar.poscar.natoms
        chgcarpath = os.path.abspath(chgcar_filename)
        data = []
        
	with open(acf_path) as f:
            print "Reading ACF"
            raw = f.readlines()
            headers = [s.lower() for s in raw.pop(0).split()]
            raw.pop(0)
        while True:
            l = raw.pop(0).strip()
            if l.startswith("-"):
                break
            vals = map(float, l.split()[1:])
            data.append(dict(zip(headers[1:], vals)))
        for l in raw:
            toks = l.strip().split(":")
            if toks[0] == "VACUUM CHARGE":
                self.vacuum_charge = float(toks[1])
            elif toks[0] == "VACUUM VOLUME":
                self.vacuum_volume = float(toks[1])
            elif toks[0] == "NUMBER OF ELECTRONS":
                self.nelectrons = float(toks[1])
        self.data = data
    
    def get_charge_transfer(self, atom_index):
        """
        Returns the charge transferred for a particular atom. Requires POTCAR
        to be supplied.
        Args:
        atom_index:
        Index of atom.
        Returns:
        Charge transfer associated with atom from the Bader analysis.
        Given by final charge on atom - nelectrons in POTCAR for
        associated atom.
        """
        if self.potcar is None:
            raise ValueError("POTCAR must be supplied in order to calculate "
                "charge transfer!")
        potcar_indices = []
        for i, v in enumerate(self.natoms):
            potcar_indices += [i] * v
        nelect = self.potcar[potcar_indices[atom_index]].nelectrons
        return self.data[atom_index]["charge"] - nelect#nelect - self.data[atom_index]["charge"] 

    
    def get_charge(self, atom_index_list= None):
        if atom_index_list:
            for i in atom_index_list:
                 print "Charge at atom index ", i, " is ", self.data[i]["charge"]
    
    def get_oxidation_state_decorated_structure(self):
        structure = self.chgcar.structure
        charges = [self.get_charge_transfer(i) for i in range(len(structure))]
        structure.add_oxidation_state_by_site(charges)
        return structure
    
class Scan_Interface(Structure):
    """
    moves the ligand
    """

    def __init__(self,structure = None, search_box= [0.1,0.1,-0.1], search_step= 0.01, energy_model='bader_coulomb', interface_barrier_height= 0.62):
	self.structure = structure
	self.search_box = search_box            #cuboid search box
        self.search_step = search_step          #step size of search
	self.energy_model = energy_model        #energy model
	self.interface_barrier_height = interface_barrier_height #dividing height for atoms to be shifted
	
     
    def ligand_shifts(self):
	vec = []
	box_a = self.search_box[0]
	box_b = self.search_box[1]
	box_c = self.search_box[2]
        step = self.search_step
	limit_a = box_a/step
	limit_b = box_b/step
	limit_c = box_c/step
	print limit_a, limit_b, limit_c
        for i in range(0, int(limit_a)):
        	for j in range(0,int(limit_b)):
            		for k in range(0,int(-1*limit_c)):
                		vec.append([i*step, j*step, k*step])

        return vec

    def get_struct_energy(self,structure):
	ce = ColoumbEnergy()
	strt = structure
    	#direct_sum_energy=ce.get_energy(strt_ox_iface) #other energy models
    	if self.energy_model == 'bader_coulomb':
		energy=ce.get_bader_coulomb_energy(strt)
    	#ewald_energy = ce.get_ewald_sum(strt_ox_iface)
	return energy 


    def trial_interface(self):
	strts = []
	energies = []
	vecs = self.ligand_shifts()
        old_struct_energy = self.get_struct_energy(self.structure)
	print "calculated old energy ... searching for lower energies"
	print old_struct_energy
	print enumerate(vecs)
	for k,j in enumerate(vecs):
        	trans = self.structure.copy()
        	for i, sitei in enumerate(trans.frac_coords):
                	if sitei[2]>self.interface_barrier_height:  #gets ligand atoms to translate sites
                		trans.translate_sites([i],j)
        	if not trans.is_valid(tol = 1.0): #tolearnce for closeness of atoms set to 1 A
                	print "possibly too close atoms"
                else:
			egy=self.get_struct_energy(structure=trans)
			print egy
                	print "energy= ", egy
                	if egy<old_struct_energy:
                       		print "Lower energy", j
				strts.append(trans)
				energies.append(egy)
				if len(strts)>10:
					return strts, energies
	return strts, energies







if __name__=='__main__':
    
#	relaxed_structure = Structure.from_file("POSCAR")
	bader_path= "./bader CHGCAR -ref Charge_sum"
	BA = Bader_Analysis(bader_path=bader_path) # optional user specifies bader path 
	relaxed_oxidated_structure = BA.get_oxidation_state_decorated_structure()
	print "got oxidated structure.. starting optimization"
	print relaxed_oxidated_structure
	lig_part_charge = []
	slab_part_charges = []
	total = []
	for i,sitei in enumerate(relaxed_oxidated_structure.as_dict()['sites']):
		total.append(sitei['species'][0]['oxidation_state'])
		if i in [0, 1, 2, 3, 4, 5, 54, 55, 56, 57, 58, 107, 108, 109,110]:
			print sitei['species'][0]['element'], sitei['species'][0]['oxidation_state']
			lig_part_charge.append(sitei['species'][0]['oxidation_state']) 
		else:
			slab_part_charges.append(sitei['species'][0]['oxidation_state'])
	#calculate total charge transfers from bader
	print "total sum: ", np.sum(total)
	print "ligand sum: ", np.sum(lig_part_charge)
	print "slab sum: ", np.sum(slab_part_charges)			
	SI = Scan_Interface(relaxed_oxidated_structure)
	trial_structs, energies = SI.trial_interface()
	dict_structs = dict(zip(energies, trial_structs))
        #trial sorting --
	print sorted(dict_structs)
	#print trial_structs, energies
	print trial_structs[0]
	print energies[0]

















#    #STEP 1 : Read Chargecar files from initial relaxation run **provide appropriate path to working directory

#    print "Reading begins.."
#    Vol_obj_1 = Chgcar.from_file("../AECCAR0")
#    print "Read AECCAR0"
#    Vol_obj_2 = Chgcar.from_file("../AECCAR2")
#    print "Read AECCAR2"


#   Chgcar_file = Chgcar.from_file("../CHGCAR")
#   print "Read CHGCAR"



#   Vol_obj_sum = Vol_obj_1.linear_add(Vol_obj_2)
#   print "Added to charge sum writing charge sum"

#   Vol_obj_sum.write_file("./Charge_sum")

    #STEP 2: Bader executable path ..run Bader 
#   exe = "~/Softwares/Bader/bader CHGCAR -ref Charge_sum"#"~/Softwares/Bader/bader ../CHGCAR -ref Charge_sum"

#   os.system(exe)

    #STEP 3: Path to Bader ACF.dat
#   path1 = "/home/josh/Research/Work/Direct_Solvation_Model/data/Raw_backup/PbS_Ligands_semiauto/PbS111_HOPO/Interface/Static/Bader_analysis/"
#   print "Reading interface"
#   BA_1 = Bader_Analysis(acf_path=path1+"/ACF.dat", chgcar_filename=path1+"../CHGCAR", potcar_filename= path1+"../POTCAR")
    #111 R Slab
#   print "decorated iface"
#   strt_ox_iface = BA_1.get_oxidation_state_decorated_structure()
#   print strt_ox_iface
#   trans = strt_ox_iface.copy()
    #compute coulomb sums

#   st = []
   # for i in strt_ox_iface.as_dict()['sites']:
   #     print i['species'][0]['element']
   #     st.append(i['species'][0]['oxidation_state'])
   # print np.sum(st)
#   ce = ColoumbEnergy()
#   direct_sum_energy=ce.get_energy(strt_ox_iface)
#   bader_sum_energy=ce.get_bader_coulomb_energy(strt_ox_iface)
#   ewald_energy = ce.get_ewald_sum(strt_ox_iface)
#   print "ewald_energy: ", ewald_energy
#   print "simple coulombic energy: ", direct_sum_energy
#   print "bader coulombic energy: ", bader_sum_energy
#   print "finished"
    #STEP 4: Scan ligand over interface for regions of lower potential [limits of x,y,z with step size]
#   vec_list = []
#   vec = []
#   step = 0.01
#   for i in range(0, 10):
#       for j in range(0,10):
#           for k in range(0,10):
#               vec.append([i*step, j*step, -1*k*step])

#   print strt_ox_iface.sites[0]
#   dcoul = []
#   bcoul = []
#   for k,j in enumerate(vec):
#       trans = strt_ox_iface.copy()
#    	for i, sitei in enumerate(trans.frac_coords):
#		if sitei[2]>0.65:  #gets ligand atoms to translate sites
#			trans.translate_sites([i],j)
	#limit of atomic distances 		
# 	if not trans.is_valid(tol = 1.0):
#		print "possibly too close atoms"
#	else:
#		filename = "./PbS111_HOPO/POSCAR"+str(k)+".vasp"
#		print filename, str(j)
#		ewd=ce.get_ewald_sum(trans) 
#		bcoulomb=ce.get_bader_coulomb_energy(trans)
#		print "ewald sum= ", ewd
#		print "bader coulomb sum= ", bcoulomb
#		if ewd<ewald_energy:
#			print "Lower energy"
#			pos = Poscar(trans)
#			pos.write_file(filename)


###########	
#	dcoul.append(ewd)

#	bcoul.append(bcoulomb)
#    plt.plot(dcoul, label='D-Coulomb')
#    plt.plot(bcoul, label='B-Coulomb')
    #plt.plot(ewd, label='Ewald')
 #   plt.savefig('Coulomb_sums.png')
	

# In[8]:




# In[ ]:



