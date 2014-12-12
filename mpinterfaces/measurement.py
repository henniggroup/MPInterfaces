"""

combines instrument, calibrate and interfaces to perform the calibration and run the actual jobs

"""
from mpinterfaces.calibrate import Calibrate

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
