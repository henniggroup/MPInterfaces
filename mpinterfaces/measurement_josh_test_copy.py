"""

combines instrument, calibrate and interfaces to perform the calibration and run the actual jobs

"""
from mpinterfaces.calibrate import Calibrate
from pymatgen.io.vasp_io import Vasprun
import sys

class Measurement(Calibrate):
    
    """
    
    Perfor the actual measurement using the optimum knob settings
    takes calibratemolecule, calibrateslab and interface objects and
    perform the interface calculations and compute binding energy, bandstructure etc. 
    
    """

    def __init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./Measurement'):
        Calibrate.__init__(self, incar, poscar, potcar, kpoints, setup_dir=setup_dir, parent_job_dir=parent_job_dir)
        self.interface = interface
        self.encut = None
        self.kpoints = None
        self.vac_spacing = None
        self.slab_thickness = None
        self.jobs = []
        self.handlers = []


    def knob_settings(self):
        """
        go through the parent dir and get all encut, kpoints and energies
        also vac spacing and slab thinckness for slab calulations
        these values willl be used to do the actual interface measurements
        use Vasprun class to get the afore mentioned values from the xml files
        should not proceed if the calculations are not done
        should update the incar, poscar, potcar, kpoints objects according to
        the knob_settings
        
        
        """

	#first begin from parent directory: . job_dir is the directory of all the ENCUT runs 

	file = open("Calibration_data", 'a')
	file.write("ENCUT or KPOINT or vacuum space in A, Energy per unit volume \n ENCUT convergence \n")
	#change to directory in sequence of folders 
	for dir_name in range(500, 800, 50): #sequence of job_dirs or entire parent directory
    		os.chdir(#directory containing ultimate vaspruns to parse)
		Delta_E = 0
		flag = 0
    		try:   
        		vasp_out = Vasprun("vasprun.xml")
        		Energy_per_vol = vasp_out.final_energy/vasp_out.lattice.volume
        		os.chdir(#to parent dir"../../")
        		file.write(str(dir_name)+', '+str(Energy_per_vol)+'\n')
        		Delta_E = Energy_per_vol - Delta_E
			if Delta_E<0.01:
				Knob_ENCUT = dir_name
				flag = 1
				self.incar["ENCUT"] = Knob_ENCUT
        		#print "writing"
    		except:
        		print "job not done, vasprun not ready to parse" 
	if flag == 0:
		file.write("ENCUT convergence not reached")
	file.write("\n KPOINTS Convergence")
	for dir_name in #job_dir seq for kpoints:
		os.chdir(#to ultimate vasprun bearing directory "Kpoints numbered")
		Delta_E = 0
		flag = 0 
		try:
			vasp_out = Vasprun("vasprun.xml")
                	Energy_per_vol = vasp_out.final_energy/vasp_out.lattice.volume
                	os.chdir(#to parent dir to write to file"../../")
                	file.write(str(dir_name)+', '+str(Energy_per_vol)+'\n')
                	Delta_E = Energy_per_vol - Delta_E
                	if Delta_E<0.01:
                		Knob_KPOINT = dir_name
                        	flag = 1
				self.kpoints = #Kpoint constructor with Knob_KPOINT as input
		except:
			print "job not done, vasprun not ready to parse"

	if flag == 0: 
		file.write("Kpoint Convergence not reached")
	pass
	for dir_name in #job_dir seq for slab vacuum:
                os.chdir(#ulimate directory containing vasprun.xml"vacuum space numbered") #path according to sequence
                Delta_E = 0
                flag = 0
                vasp_out = Vasprun("vasprun.xml")
                Energy_per_vol = vasp_out.final_energy/vasp_out.lattice.volume
                os.chdir(#to parent directory"../../")
                file.write(str(dir_name)+', '+str(Energy_per_vol)+'\n')
                Delta_E = Energy_per_vol - Delta_E
                if Delta_E<0.01:
                        Knob_Vacuum = dir_name
                        flag = 1
			self.vac_spacing = #
                	self.min_thickness =#
                	#construct slab with vac_spacing and min_thickness, write to structure
                	self.poscar =#Poscar from structure	
        if flag == 0:
                file.write("Vacuum spacing converence Convergence not reached")
  
                pass

    def setup(self, interfaces):
        """
        create a list of jobs to run
        exampls:- change coverage, site occupancies,
        interfaces: list of interfaces for which jobs must be created
        also run molecule, slab relaxation calcs
        """
        for iface in interfaces:
            job_dir  = self.parent_job_dir + iface.name
            self.job_dirs.append(job_dir)
            vis = MPINTVaspInputSet(iface.name, self.incar, self.poscar, self.potcar, self.kpoints)
            #the job command can be overrridden in the run method
            job = MPINTVaspJob(["pwd"], final = True, setup_dir=self.setup_dir, job_dir=job_dir, vis=vis, auto_npar=False, auto_gamma=False)
            self.jobs.append(job)

    def make_measurements(self):
        """
        To calculate binding energy of interface.
        Needs to get relaxed molecule energy, relaxed slab energy and relaxed interface energies
        Get the molecule and slab energies from a relaxation run
        of optimum paramter set for Molecule and Slab (when is
        that relaxation run done? NOTE: static done with
        as of Calibration, do before this call)
        Use Vasprun() to get energies, same method as knob_settings in Measurement.

        Goes to the required directories for Molecule and Slab and uses Vasprun to parse through the vasprun.xml and gets the energy, similar for 
        """
        
        pass
        



#test
if __name__=='__main__':
#    interfaces = Interface objects    
#    measure = Measurement(incar, poscar, potcar, kpoints)
#    measure.knob_settings()
#    measure.setup(interfaces)
#    measure.run(['qsub','job_script'])
