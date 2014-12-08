"""
The instrument module

"""


import sys
import os, shutil
import shlex, subprocess
import time
import datetime
from pprint import pprint
import logging

import numpy as np

from pymatgen import Composition, Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.design_patterns import Enum
#from pymatgen.io.smartio import read_structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vaspio_set import DictVaspInputSet #MPGGAVaspInputSet

#from custodian.vasp.handlers import VaspErrorHandler, FrozenJobErrorHandler, MeshSymmetryErrorHandler, NonConvergingErrorHandler
from custodian.custodian import Custodian, Job, gzip_dir
from custodian.vasp.jobs import VaspJob
from custodian.ansible.actions import FileActions, DictActions
from custodian.vasp.interpreter import VaspModder


VASP_INPUT_FILES = {"INCAR", "POSCAR", "POTCAR", "KPOINTS"}

VASP_OUTPUT_FILES = ['DOSCAR', 'INCAR', 'KPOINTS', 'POSCAR', 'PROCAR',
                     'vasprun.xml', 'CHGCAR', 'CHG', 'EIGENVAL', 'OSZICAR',
                     'WAVECAR', 'CONTCAR', 'IBZKPT', 'OUTCAR']


class myVaspInputSet(DictVaspInputSet):
    """
    subcalss DictVaspInputSet and customize write_input method
    myVIS.yaml should be in the MODULE_DIR(if thats being used)
    use user_incar_settings to override the defaults in myVIS.yaml
    
    """
    def __init__(self, name, incar, poscar, potcar, kpoints,**kwargs ):#config_dict, user_incar_settings=None, **kwargs):
        """
        default INCAR from config_dict
        
        """
        self.name = name
        
        self.incar = Incar.from_dict(incar.as_dict())#copy.deepcopy(incar)
        self.poscar = Poscar.from_dict(poscar.as_dict())
        self.potcar = Potcar.from_dict(potcar.as_dict())
        self.kpoints = Kpoints.from_dict(kpoints.as_dict())
        
        config_dict = {}
        config_dict['INCAR'] = self.incar.as_dict()
        config_dict['POTCAR'] = dict(zip(self.potcar.as_dict()['symbols'], self.potcar.as_dict()['symbols'])) #caution the key and the value are not always the same
        config_dict['KPOINTS'] = self.kpoints #kpoints.as_dict()
        #self.user_incar_settings = self.incar.as_dict()        
        
        DictVaspInputSet.__init__(self, name, config_dict, ediff_per_atom=False, **kwargs)

    def write_input(self, job_dir, make_dir_if_not_present=True, write_cif=False):
        """
        the input files are written to the job_dir
        process(if needed) and write the input files in each directory
        structures read from the poscar files in the directory
        
        """
        d = './'+job_dir #self.name        
        if make_dir_if_not_present and not os.path.exists(d):
            os.makedirs(d)
        self.incar.write_file(os.path.join(d, 'INCAR'))
        self.kpoints.write_file(os.path.join(d, 'KPOINTS'))
        self.potcar.write_file(os.path.join(d, 'POTCAR'))
        self.poscar.write_file(os.path.join(d, 'POSCAR'))
        #print os.getcwd()

    def as_dict(self):
        d = super(myVaspInputSet, self).as_dict()
        return d



class myVaspJob(Job):
    """
    customize the VASPJob class: setup, run and postprocess functions overridden
    
    Args:
       vasp_cmd : the command to be issued in wach job_dir
       setup_dir : directory that has the setup files for creating the rest of the vasp inputs
       job_dir : the directory from which the jobs will be launched
    
    """
    def __init__(self, vasp_cmd, output_file="vasp.out", setup_dir='.', job_dir='untitled', suffix="",
                 final=True, gzipped=False, backup=True,
                 vis=None, auto_npar=True,
                 auto_gamma=True, settings_override=None,
                 gamma_vasp_cmd=None, copy_magmom=False):

        self.vasp_cmd = vasp_cmd
        self.output_file = output_file
        self.setup_dir = setup_dir
        self.job_dir = job_dir
        self.final = final
        self.backup = backup
        self.gzipped = gzipped
        self.vis = vis
        self.suffix = suffix
        self.settings_override = settings_override
        self.auto_npar = auto_npar
        self.auto_gamma = auto_gamma
        self.gamma_vasp_cmd = gamma_vasp_cmd
        self.copy_magmom = copy_magmom

    def setup(self):
        """
        looks for the set up files(POSCAR, submit_job etc) in the setup_dir
        uses those files to create the vasp input set in the job_dir
        the current setup looks only for the poscar file in the setup directory
        """
        files = os.listdir(self.setup_dir)
        num_structures = 0
        if not set(files).issuperset(VASP_INPUT_FILES):
            for f in files:
                try:
                    struct = Structure.from_file(f)
                    num_structures += 1
                except:
                    pass
            if num_structures != 1:
                raise RuntimeError("{} structures found. Unable to continue."
                                   .format(num_structures))
            else:
                self.vis.write_input(self.job_dir)
                #shutil.copy('submit_job', self.job_dir)

        if self.backup:
            os.chdir(self.job_dir)
            for f in VASP_INPUT_FILES:
                shutil.copy(f, "{}.orig".format(f))
            os.chdir('../')                

        #if self.settings_override is not None:
        #    VaspModder().apply_actions(self.settings_override)
        #os.chdir('../')            

    def run(self):
        """
        move to the job_dir, launch the job and get out
        """
        os.chdir(self.job_dir)
        cmd = list(self.vasp_cmd)
        if self.auto_gamma:
            vi = VaspInput.from_directory(".")
            kpts = vi["KPOINTS"]
            if kpts.style == "Gamma" and tuple(kpts.kpts[0]) == (1, 1, 1):
                if self.gamma_vasp_cmd is not None and which(
                        self.gamma_vasp_cmd[-1]):
                    cmd = self.gamma_vasp_cmd
                elif which(cmd[-1] + ".gamma"):
                    cmd[-1] += ".gamma"
        logging.info("Running {}".format(" ".join(cmd)))
        with open(self.output_file, 'w') as f:
            p = subprocess.Popen(cmd, stdout=f)
            os.chdir('../')
        return p

    def postprocess(self):
        ## major clean up
        for f in VASP_OUTPUT_FILES + [self.output_file]:
            if os.path.exists(f):
                if self.final and self.suffix != "":
                    shutil.move(f, "{}{}".format(f, self.suffix))
                elif self.suffix != "":
                    shutil.copy(f, "{}{}".format(f, self.suffix))

        if self.copy_magmom and not self.final:
            try:
                outcar = Outcar("OUTCAR")
                magmom = [m['tot'] for m in outcar.magnetization]
                incar = Incar.from_file("INCAR")
                incar['MAGMOM'] = magmom
                incar.write_file("INCAR")
            except:
                logging.error('MAGMOM copy from OUTCAR to INCAR failed')

        if self.gzipped:
            gzip_dir(".")

    def name(self):
         return self.__class__.__name__
    


class Calibrate(object):
    
    """
        
    The base class for creating vasp work flows for calibrating the input parameters for different systems
    
    
    """

    def __init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./'):
        self.parent_job_dir = parent_job_dir
        self.setup_dir = setup_dir
        self.incar = incar
        self.poscar = poscar
        self.potcar =potcar
        self.kpoints = kpoints
        self.vis = []
        self.jobs = []
        self.handlers = [ ] #example:- handlers = [VaspErrorHandler(), FrozenJobErrorHandler(), MeshSymmetryErrorHandler(), NonConvergingErrorHandler()]

    def encut_cnvg(self, encut_list):
        """
        create ENCUT convergence workflow
        """
        for encut in encut_list:
            print 'ENCUT = ', encut
            self.job_dir = self.parent_job_dir + str(encut)
            self.incar['ENCUT'] = encut
            vis = myVaspInputSet('encut_'+str(encut), self.incar, self.poscar, self.potcar, self.kpoints)
            job = myVaspJob(["pwd"], final = True, setup_dir=self.setup_dir, job_dir=self.job_dir, vis=vis, auto_npar=False, auto_gamma=False)
            self.jobs.append(job)

    def kpoints_cnvg(self, Grid_type = 'M', kpoints_tuple = ((7, 7, 7), (10, 10, 10)), conv_step = 1):
        """
        must be implemented in the deriving class
        """
        pass

        
    def run(self):
        """
        run the vasp jobs through custodian
        set up the custodian that we want to run : just a single run
        just one vasp job that uses the above defined input set
        first param in VaspJob = vasp_cmd eg:- ["mpirun", "pvasp.5.2.11"]
        create multiple vaspjob object, each with differnt dictvaspinputset object
        consider subclassing VaspJob and overriding the post_process method
        create a custodian task using the jobs and error handlers
        consider subclassing Custodian and overriding _run_job method to modify when _do_check method is called
        it must be called when the calcualtion is done for example: check the OUTCAR and check for
        'writing wavefunctions' or maybe not

        disable error handlers on the job and the postprocess since the 'job' here just submitting the job to the queue 

        """
        c_params = {'jobs': [j.as_dict() for j in self.jobs], 'handlers': [h.as_dict() for h in self.handlers], 'max_errors': 5}
        c = Custodian(self.handlers, self.jobs, max_errors=5)
        c.run()


        
class CalibrateMolecule(Calibrate):
    
    """
    
    Calibrate paramters for Molecule calculations
    
    """
    
    def __init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./'):
        Calibrate.__init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./')

    def kpoints_cnvg(self, Grid_type = 'M', kpoints_tuple = ((7, 7, 7), (10, 10, 10)), conv_step = 1):
        """
        If molecule constructs single kpoint file with MP 1x1x1.
        Calls Kpoints constructors according to Grid_type: Monkhorst Pack Automatic is default,
        G for Gamma centered automatic; ONLY MP method implemented now others can be added 
        if_slab flag is for the application of constraint to z-axis for slab,
        kpoints_list describes the start and end of kpoints set
        eg: user can pass ((6, 6, 6), ((10, 10, 10))) and conv_step for a convergence to be done for 6x6x6 to 10x10x10,
        if_slab can be switched to True constraints z to default to 1, that is 6x6x1 to 10x10x1
        first defines the list of kpoints according to the user input,
        user needs to give only the start and end kpoint, whether it is for a slab, molecule or simple bulk
        """
        kpoint = (1, 1, 1)
        self.kpoints = self.kpoints.monkhorst_automatic(kpts=kpoint)
        K = list(kpoint)
        print 'KPOINTmesh = ', str(K[0])+'x'+str(K[1])+'x'+str(K[2])
        self.job_dir = self.parent_job_dir + str(K[0])+'x'+str(K[1])+'x'+str(K[2])
        vis = myVaspInputSet('kpoint_'+str(K[0])+'x'+str(K[1])+'x'+str(K[2]), self.incar, self.poscar, self.potcar, self.kpoints)
        job = myVaspJob(["pwd"], final = True, setup_dir=self.setup_dir, job_dir=self.job_dir, vis=vis, auto_npar=False, auto_gamma=False)
        self.jobs.append(job)
        #return "kpoints constructed for simple molecule"
        
                

class CalibrateBulk(Calibrate):
    
    """
    
    Calibrate paramters for Bulk calculations
    
    """
    
    def __init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./'):
        Calibrate.__init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./')

        
    def kpoints_cnvg(self, Grid_type = 'M', kpoints_tuple = ((7, 7, 7), (10, 10, 10)), conv_step = 1):
        if Grid_type == 'M':
            print "at Grid"
            conv_list = list(kpoints_tuple) #local list convergence_list , convert from tuple because constructor takes tuple as argument
            start = list(conv_list[0])
            end = list(conv_list[1])
            if (conv_step):
                for x in range(1+start[0], end[0], conv_step):
                    conv_list.append([x, x, x])                
                for kpoint in conv_list:
                    self.kpoints = self.kpoints.monkhorst_automatic(kpts = kpoint)
                    K = list(kpoint)
                    print 'KPOINTmesh = ', str(K[0])+'x'+str(K[1])+'x'+str(K[2])
                    self.job_dir = self.parent_job_dir + str(K[0])+'x'+str(K[1])+'x'+str(K[2])
                    vis = myVaspInputSet('kpoint_'+str(K[0])+'x'+str(K[1])+'x'+str(K[2]), self.incar, self.poscar, self.potcar, self.kpoints)
                    job = myVaspJob(["pwd"], final = True, setup_dir=self.setup_dir, job_dir=self.job_dir, vis=vis, auto_npar=False, auto_gamma=False)
                    self.jobs.append(job)
        

                    
class CalibrateSlab(Calibrate):
    
    """
    
    Calibrate paramters for Slab calculations
    
    """
    
    def __init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./'):
        Calibrate.__init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./')

        
    def kpoints_cnvg(self, Grid_type = 'M', kpoints_tuple = ((7, 7, 7), (10, 10, 10)), conv_step = 1):
        if Grid_type == 'M':
            conv_list = list(kpoints_tuple) #local list convergence_list , convert from tuple because constructor takes tuple as argument
            start = list(conv_list[0])
            end = list(conv_list[1])
            if start[2] != 1 or end[2] != 1:
                return "Kpoints not for slab input!"
            elif (conv_step):
                for x in range(1+start[0], end[0], conv_step):
                    conv_list.append([x, x, 1])
                for kpoint in conv_list:
                    self.kpoints = self.kpoints.monkhorst_automatic(kpts = kpoint)
                    K = list(kpoint)
                    print 'KPOINTmesh = ', str(K[0])+'x'+str(K[1])+'x'+str(K[2])
                    self.job_dir = self.parent_job_dir + str(K[0])+'x'+str(K[1])+'x'+str(K[2])
                    vis = myVaspInputSet('kpoint_'+str(K[0])+'x'+str(K[1])+'x'+str(K[2]), self.incar, self.poscar, self.potcar, self.kpoints)
                    job = myVaspJob(["pwd"], final = True, setup_dir=self.setup_dir, job_dir=self.job_dir, vis=vis, auto_npar=False, auto_gamma=False)
                    self.jobs.append(job)
	            

                    
class CalibrateInterface(Calibrate):
    
    """
    
    Calibrate paramters for Interface(= Slab+ligand+solvent) calculations
    
    """    
    
    def __init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./'):
        Calibrate.__init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./')



class Knobs(object):
    
    """
    
    Fetch the optimum param values
    
    """
    pass



class Experiment(object):
    
    """
    
    Run the actual experiment using the optimum knob settings
    example: slab+ligand calculations
    
    """
    pass



#test
if __name__ == '__main__':

    system = 'Pt bulk'

    atoms = ['Pt']
    
    a0 = 3.965
    lvec = [ [0.5, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5] ]
    lvec = np.array(lvec) * a0
    lattice = Lattice(lvec)#.from_parameters(3.866, 3.866, 3.866, 60, 60, 60)
    structure = Structure( lattice, atoms, [ [0.0, 0.0, 0.0] ], coords_are_cartesian=False, site_properties={"magmom":[0]} )
    #a POSCAR file should be in the setup directory run directory
    structure.to(fmt="poscar", filename="POSCAR")

    incarparams = {'System':'test', 'ENCUT': 400, 'ISMEAR': 1, 'SIGMA': 0.1, 'EDIFF':1E-6}
    incar = Incar(params=incarparams)
    poscar = Poscar(structure, comment=system, selective_dynamics=None, true_names=True, velocities=None, predictor_corrector=None)
    atoms = poscar.site_symbols
    potcar = Potcar(symbols=atoms, functional='PBE', sym_potcar_map=None)
    kpoints = Kpoints.monkhorst_automatic(kpts=(16, 16, 16), shift=(0, 0, 0))#{'grid_density': 1000} #

#    calmol = CalibrateMolecule(incar, poscar, potcar, kpoints)
    calbulk = CalibrateBulk(incar, poscar, potcar, kpoints)    
    calbulk.encut_cnvg(range(400,800,100))
    calbulk.kpoints_cnvg(kpoints_tuple = ((7, 7, 7), (11, 11, 11))) #for bulk
    calbulk.run()



#note: write the default yaml to the directory where the jobs are run. This is useful later on for comparing different job runs in that direcctory
