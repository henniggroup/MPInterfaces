"""

combines instrument, calibrate and interfaces to perform the calibration and run the actual jobs

"""

import numpy as np
from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from mpinterfaces.calibrate import Calibrate, CalibrateMolecule, CalibrateSlab
from mpinterfaces.interface import Interface
        

class Measurement(object):
    
    """
    
    Perfor the actual measurement using the optimum knob settings
    takes calibratemolecule, and calibrateslab
    setups the interface using the poscars of calmol and calslab objects
    also uses the knob_settings of calmol and calslab to set up the input set
    for the interface calulation
    perform the interface calculations and compute binding energy, bandstructure etc. 
    
    """

    def __init__(self, cal_objs, setup_dir='.', parent_job_dir='.', job_dir='./Measurement'):
        #self.interface = interface
        self.encut = None
        self.kpoints = None
        self.vac_spacing = None
        self.slab_thickness = None
        self.jobs = []
        self.handlers = []
        self.calmol = []
        self.calslab = []
        self.calbulk = []                
        for obj in cal_objs:
            if isinstance(obj, CalibrateMolecule):
                self.calmol.append(obj)
            elif isinstance(obj, CalibrateSlab):
                self.calslab.append(obj)
            elif isinstance(obj, CalibrateBulk):
                self.calbulk.append(obj)                
            


    def setup(self):
        """
        uses the calmol and calslab objects to creat the interface poscar
        and other input files for vasp calculation
        create a list of jobs to run
        exampls:- change coverage, site occupancies,
        also run molecule, slab relaxation calcs
        """
        #job_dir  = self.parent_job_dir + iface.name
        #self.job_dirs.append(job_dir)
        #vis = MPINTVaspInputSet(iface.name, self.incar, self.poscar, self.potcar, self.kpoints)
            #the job command can be overrridden in the run method
        #job = MPINTVaspJob(["pwd"], final = True, setup_dir=self.setup_dir, job_dir=job_dir, vis=vis, auto_npar=False, auto_gamma=False)
        #self.jobs.append(job)
        pass

    def run(self, job_cmd=None):
        #override the job_cmd if provided
        if job_cmd :
            for j in self.jobs:
                j.job_cmd = job_cmd
                
        c_params = {'jobs': [j.as_dict() for j in self.jobs], 'handlers': [h.as_dict() for h in self.handlers], 'max_errors': 5}
        c = Custodian(self.handlers, self.jobs, max_errors=5)
        c.run()
        

    def make_measurements(self):
        """
        To calculate binding energy of interface.
        Needs to get relaxed molecule energy, relaxed slab energy and relaxed interface energies
        Get the molecule and slab energies from a relaxation run
        of optimum paramter set for Molecule and Slab.
        Goes to the required directories for Molecule and Slab and uses Vasprun to parse through the vasprun.xml and gets the energy,
        and compute the required quantities
        """
        
        pass



#test
if __name__=='__main__':
    system = 'Pt bulk'
    atoms = ['Pt']
    a0 = 3.965
    lvec = [ [0.5, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5] ]
    lvec = np.array(lvec) * a0
    lattice = Lattice(lvec)
    structure = Structure( lattice, atoms, [ [0.0, 0.0, 0.0] ],
                           coords_are_cartesian=False,
                           site_properties={"magmom":[0]} )

    incarparams = {'System':'test',
                   'ENCUT': 400,
                   'ISMEAR': 1,
                   'SIGMA': 0.1,
                   'EDIFF':1E-6}
    incar = Incar(params=incarparams)
    poscar = Poscar(structure, comment=system,
                    selective_dynamics=None,
                    true_names=True, velocities=None, predictor_corrector=None)
    potcar = Potcar(symbols=poscar.site_symbols, functional='PBE', sym_potcar_map=None)
    kpoints = Kpoints.monkhorst_automatic(kpts=(16, 16, 16), shift=(0, 0, 0))

#    measure = Measurement([calmol, calslab], job_dir='./Measurement')

