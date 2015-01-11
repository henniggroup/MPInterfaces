"""

combines instrument, calibrate and interfaces to perform the calibration
and run the actual jobs

"""
import shutil
import os
import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints

from mpinterfaces.calibrate import Calibrate, CalibrateMolecule,\
      CalibrateSlab, CalibrateBulk
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
    def __init__(self, cal_objs, setup_dir='.', parent_job_dir='.',
                 job_dir='./Measurement'):
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
        self.cal_objs = cal_objs
        self.job_dir = job_dir
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
        #vis = MPINTVaspInputSet(iface.name, self.incar, self.poscar,
            #self.potcar, self.kpoints)
            #the job command can be overrridden in the run method
        #job = MPINTVaspJob(["pwd"], final = True, setup_dir=self.setup_dir,
            # job_dir=job_dir, vis=vis, auto_npar=False, auto_gamma=False)
        #self.jobs.append(job)
        pass

    def run(self, job_cmd=None):
        """ run jobs """
        for cal in self.cal_objs:
            if cal.calc_done:
                cal.run()
            elif not cal.isrunning:
                cal.setup()
                cal.run()
            else:
                print 'calibration calc still running'
                print 'try again later'

    def setup_static_job(self, cal):
        """
        setup static jobs for the calibrate objects
        copies CONTCAR to POSCAR
        and
        set NSW = 0
        """
        if cal.calc_done:
            job_dir = self.job_dir+os.sep+'STATIC'
            contcar_file = cal.parent_job_dir+os.sep+cal.job_dir+os.sep+'CONTCAR'            
            cal.poscar = Poscar.from_file(contcar_file)
            cal.incar['NSW'] = 0
            cal.add_job(job_dir=job_dir)
        else:
            cal.jobs = []
            print 'previous calc in the dir, ', cal.job_dir, 'not done yet or is still running'
            print 'Not setting up the measurement job\n'            

    def setup_solvation_job(self, cal):
        """
        setup solvation jobs for the calibrate objects
        copies WAVECAR
        and
        sets the solvation params in the incar file
        """
        if cal.calc_done:       
            job_dir = self.job_dir+os.sep+'SOL'       
            cal.incar['LSOL'] = '.TRUE.'
            cal.incar['EB_K'] = 80
            if not os.path.exists(job_dir):            
                os.makedirs(job_dir)
            wavecar_file = cal.parent_job_dir+os.sep+cal.job_dir+os.sep+'WAVECAR'
            shutil.copy(wavecar_file, job_dir+os.sep+'WAVECAR')
            cal.add_job(job_dir=job_dir)
        else:
            cal.jobs = []
            print 'previous calc in the dir, ', cal.job_dir, 'not done yet or is still running'
            print 'Not setting up the measurement job\n' 
        
    def make_measurements(self):
        """
        To calculate binding energy of interface.
        Needs to get relaxed molecule energy,
        relaxed slab energy and relaxed interface energies
        Get the molecule and slab energies from a relaxation run
        of optimum paramter set for Molecule and Slab.
        Goes to the required directories for Molecule and Slab
        and uses Vasprun to parse through the vasprun.xml and gets the energy,
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
                    true_names=True, velocities=None,
                    predictor_corrector=None)
    potcar = Potcar(symbols=poscar.site_symbols, functional='PBE',
                    sym_potcar_map=None)
    kpoints = Kpoints.monkhorst_automatic(kpts=(16, 16, 16), shift=(0, 0, 0))

    cal = CalibrateBulk(incar, poscar, potcar, kpoints,
                        job_dir='test', job_cmd=['ls','-lt'])
    #list of calibrate objects
    cal_objs = [cal]
    #check whether the cal jobs were done 
    Calibrate.check_calcs(cal_objs)
    #set the measurement
    measure = Measurement(cal_objs, job_dir='./Measurements')
    #set the measurement jobs
    for cal in cal_objs:                    
        measure.setup_static_job(cal)
        measure.setup_solvation_job(cal)
    measure.run()

