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
from pymatgen.io.smartio import read_structure
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


#---------------------------------------------------------------------------------------------------


#subcalss DictVaspInputSet and customize write_input method
#myVIS.yaml should be in the MODULE_DIR
#use user_incar_settings to override the defaults in myVIS.yaml
class myVaspInputSet(DictVaspInputSet):
    
    def __init__(self, name, config_dict, user_incar_settings=None, **kwargs):
        #default INCAR from config_dict
        defaults = {'IMAGES': 1, 'IBRION': 1, 'NFREE': 2, 'ISYM': 0,
                    'LORBIT': 0, 'LCHARG': False}
        if user_incar_settings:
            defaults.update(user_incar_settings)
        DictVaspInputSet.__init__(self, name, config_dict, user_incar_settings=defaults, ediff_per_atom=False, **kwargs)
        self.name = name

    #default output_dir is the current directory, directory for each vasp input set can be set here
    #in this example self.name is used to name the directories
    #process(if needed) and write the input files in each directory
    #structures read from the poscar files in the directory
    def write_input(self, structures, output_dir, make_dir_if_not_present=True, write_cif=False):
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        d = '.'#+self.name
        s = structures
        self.get_incar(s).write_file(os.path.join(d, 'INCAR'))
        self.get_kpoints(s).write_file(os.path.join(d, 'KPOINTS'))
        self.get_potcar(s).write_file(os.path.join(d, 'POTCAR'))
        self.get_poscar(s).write_file(os.path.join(d, 'POSCAR'))

#---------------------------------------------------------------------------------------------------

#customize the VASPJob class: setup, run and postprocess functions overridden
class myVaspJob(VaspJob):
    def __init__(self, vasp_cmd, output_file="vasp.out", job_dir='.', suffix="",
                 final=True, gzipped=False, backup=True,
                 default_vasp_input_set=None, auto_npar=True,
                 auto_gamma=True, settings_override=None,
                 gamma_vasp_cmd=None, copy_magmom=False):

        self.vasp_cmd = vasp_cmd
        self.output_file = output_file
        self.job_dir = job_dir
        self.final = final
        self.backup = backup
        self.gzipped = gzipped
        self.default_vis = default_vasp_input_set
        self.suffix = suffix
        self.settings_override = settings_override
        self.auto_npar = auto_npar
        self.auto_gamma = auto_gamma
        self.gamma_vasp_cmd = gamma_vasp_cmd
        self.copy_magmom = copy_magmom

    def setup(self):
        files = os.listdir(".")
        #print files
        num_structures = 0
        if not set(files).issuperset(VASP_INPUT_FILES):
            for f in files:
                try:
                    struct = read_structure(f)
                    num_structures += 1
                except:
                    pass
            if num_structures != 1:
                raise RuntimeError("{} structures found. Unable to continue."
                                   .format(num_structures))
            else:
                os.mkdir(self.job_dir)
                os.chdir(self.job_dir)                
                self.default_vis.write_input(struct, '.')

        if self.backup:
            for f in VASP_INPUT_FILES:
                shutil.copy(f, "{}.orig".format(f))

        if self.auto_npar:
            try:
                incar = Incar.from_file("INCAR")
                #Only optimized NPAR for non-HF and non-RPA calculations.
                if not (incar.get("LHFCALC") or incar.get("LRPA") or
                        incar.get("LEPSILON")):
                    if incar.get("IBRION") in [5, 6, 7, 8]:
                        # NPAR should not be set for Hessian matrix
                        # calculations, whether in DFPT or otherwise.
                        del incar["NPAR"]
                    else:
                        import multiprocessing
                        ncores = multiprocessing.cpu_count()
                        for npar in range(int(round(math.sqrt(ncores))),
                                          ncores):
                            if ncores % npar == 0:
                                incar["NPAR"] = npar
                                break
                    incar.write_file("INCAR")
            except:
                pass

        if self.settings_override is not None:
            VaspModder().apply_actions(self.settings_override)
        os.chdir('../')            

    def run(self):
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
    
        

#---------------------------------------------------------------------------------------------------

#return the vasp input set object: defines the method for creating the input files : INCAR, POSCAR, POTCAR, KPOINTS

def get_vis(name, structure, incar, poscar, potcar, kpoints):
    #create vasp input set
    #config_dict can also be obtained from yaml file
    #config_dict = loadfn(os.path.join(MODULE_DIR,"myVIS.yaml"))
    config_dict = {}
    config_dict['INCAR'] = incar.as_dict()
    config_dict['POTCAR'] = dict(zip(potcar.as_dict()['symbols'], potcar.as_dict()['symbols']))
    config_dict['KPOINTS'] = kpoints #kpoints.as_dict()
    #if you need to override the default incar settings, do it here
    incar_override =  {'IMAGES': 2, 'IBRION': 1, 'NFREE': 2, 'ISYM': 0,
                    'LORBIT': 0, 'LCHARG': False}

    #first param = name, name for the input set, used to name the directories in which the input files will be written 

    return myVaspInputSet(name, config_dict, user_incar_settings=incar_override)

#---------------------------------------------------------------------------------------------------

#create a work flow

def create_wf(structure):
    poscar = Poscar(structure, comment=system, selective_dynamics=None, true_names=True, velocities=None, predictor_corrector=None)
    atoms = poscar.site_symbols
    potcar = Potcar(symbols=atoms, functional='PBE', sym_potcar_map=None)
    kpoints = {'grid_density': 1000} #Kpoints.monkhorst_automatic(kpts=(16, 16, 16), shift=(0, 0, 0))

    jobs = []

    #Set the error handlers : no handlers set
    #example:- handlers = [VaspErrorHandler(), FrozenJobErrorHandler(), MeshSymmetryErrorHandler(), NonConvergingErrorHandler()]
 
    handlers = []

    #create the jobs list
    #this example creates vaspinputset objects with different settings fro ENCUT in incar file. the other 3 input files remain the same
    for i in [400, 500, 600]:
        incarparams = {'System':'test', 'ENCUT': i, 'ISMEAR': 1, 'SIGMA': 0.1, 'EDIFF':1E-6}
        incar = Incar(params=incarparams)
        myvis = get_vis(str(i), structure, incar, poscar, potcar, kpoints)

        #set up the custodian that we want to run : just a single run
        #just one vasp job that uses the above defined input set
        #first param in VaspJob = vasp_cmd eg:- ["mpirun", "pvasp.5.2.11"]
        #create multiple vaspjob obkject, each with differnt dictvaspinputset object
        #consider subclassing VaspJob and overriding the post_process method

        #jobs.append( myJob(str(i), 'vasp.out') )
        jobs.append( myVaspJob(["pwd"], final = True, job_dir=str(i), backup=False, default_vasp_input_set=myvis, auto_npar=False, auto_gamma=False) )
        #jobs.append( myJob2('../', 'vasp.out') )        

    #input to custodian

    c_params = {'jobs': [j.as_dict() for j in jobs], 'handlers': [h.as_dict() for h in handlers], 'max_errors': 5}

    #create a custodian task using the jobs and error handlers
    #conside subclassing Custodian and overriding _run_job method to modify when _do_check method is called
    # it must be called when the calcualtion is done for example: check the OUTCAR and check for
    # 'writing wavefunctions' or maybe not

    c = Custodian(handlers, jobs, max_errors=5)
    
    #run each job, error handlers on the job and the postprocess

    c.run()

#---------------------------------------------------------------------------------------------------

#example: 
if __name__ == '__main__':

    system = 'Pt bulk'

    atoms = ['Pt']
    
    a0 = 3.965
    lvec = [ [0.5, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5] ]
    lvec = np.array(lvec) * a0
    lattice = Lattice(lvec)#.from_parameters(3.866, 3.866, 3.866, 60, 60, 60)
    structure = Structure( lattice, atoms, [ [0.0, 0.0, 0.0] ], coords_are_cartesian=False, site_properties={"magmom":[0]} )

    #a POSCAR file should be in the parent run directory
    structure.to(fmt="poscar", filename="POSCAR")

    #create a workflow(multiple jobs, each with different incar, potcar and/or kpoints settings) and run the jobs
    create_wf(structure)

    
#---------------------------------------------------------------------------------------------------

#note: write the default yaml to the directory where the jobs are run. This is useful later on for comparing different job runs in that direcctory
