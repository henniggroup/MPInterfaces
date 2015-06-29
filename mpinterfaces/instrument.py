from __future__ import division, unicode_literals, print_function

"""
The instrument module:
defines the inputset and the job

"""

import sys
import os, shutil
import subprocess
import logging

import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vaspio_set import DictVaspInputSet #MPGGAVaspInputSet

from custodian.custodian import Job, gzip_dir, ErrorHandler
from custodian.vasp.interpreter import VaspModder

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


class MPINTVaspInputSet(DictVaspInputSet):
    """
    defines the set of input required for a vasp job i.e
    create INCAR, POSCAR, POTCAR & KPOINTS files
    
    subcalss DictVaspInputSet and customize write_input method
    myVIS.yaml should be in the MODULE_DIR(if thats being used)
    use user_incar_settings to override the defaults in myVIS.yaml
    
    """
    def __init__(self, name, incar, poscar, potcar, kpoints,
                 qadapter=None, **kwargs ):
        #config_dict, user_incar_settings=None, **kwargs):
        """
        default INCAR from config_dict
        
        """
        self.name = name
        self.incar = Incar.from_dict(incar.as_dict())
        self.poscar = Poscar.from_dict(poscar.as_dict())
        self.potcar = Potcar.from_dict(potcar.as_dict())
        self.kpoints = Kpoints.from_dict(kpoints.as_dict())
        if qadapter is not None:
            self.qadapter = qadapter.from_dict(qadapter.to_dict())
        else:
            self.qadapter = None
        
        config_dict = {}
        config_dict['INCAR'] = self.incar.as_dict()
         #caution the key and the value are not always the same        
        config_dict['POTCAR'] = dict(zip(self.potcar.as_dict()['symbols'],
                                         self.potcar.as_dict()['symbols']))
        config_dict['KPOINTS'] = self.kpoints #kpoints.as_dict()
        #self.user_incar_settings = self.incar.as_dict()        
        
        DictVaspInputSet.__init__(self, name, config_dict,
                                   ediff_per_atom=False, **kwargs)

        
    def write_input(self, job_dir, make_dir_if_not_present=True,
                     write_cif=False):
        """
        the input files are written to the job_dir
        process(if needed) and write the input files in each directory
        structures read from the poscar files in the directory
        
        """
        d = job_dir
        if make_dir_if_not_present and not os.path.exists(d):
            os.makedirs(d)
        logger.info('writing inputset to : '+d)
        self.incar.write_file(os.path.join(d, 'INCAR'))
        self.kpoints.write_file(os.path.join(d, 'KPOINTS'))
        self.potcar.write_file(os.path.join(d, 'POTCAR'))
        self.poscar.write_file(os.path.join(d, 'POSCAR'))
        if self.qadapter is not None:
            self.script_name = 'submit_script'
            with open(os.path.join(d, self.script_name), 'w') as f:
                queue_script = self.qadapter.get_script_str(job_dir)
                f.write(queue_script)           
            
    def as_dict(self):
        d = super(MPINTVaspInputSet, self).as_dict()
        return d

    
class MPINTVaspJob(Job):
    """
    defines a vasp job i.e setup the required input files and lanuch the job
    
    Args:
       job_cmd : a list, the command to be issued in each job_dir
       eg: ['qsub', 'submit_job']
       setup_dir : directory that has the setup files for creating the
       rest of the vasp inputs
       job_dir : the directory from which the jobs will be launched
    """
    def __init__(self, job_cmd, name='noname',output_file="job.out", setup_dir='.',
                 parent_job_dir='.', job_dir='untitled', suffix="",
                 final=True, gzipped=False, backup=False,
                 vis=None, auto_npar=True,
                 auto_gamma=True, settings_override=None,
                 gamma_vasp_cmd=None, copy_magmom=False, wait=True):

        self.job_cmd = job_cmd
        self.name = name
        self.output_file = output_file
        self.setup_dir = setup_dir
        self.parent_job_dir = parent_job_dir
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
        self.wait = wait

    def setup(self):
        """
        looks for the set up files(POSCAR, submit_job etc) in the setup_dir
        uses those files to create the vasp input set in the job_dir
        the current setup looks only for the poscar file in the setup directory
        """
        self.vis.write_input(self.job_dir)
        if self.backup:
            os.chdir(os.path.abspath(self.job_dir))
            for f in os.listdir('.'):
                shutil.copy(f, "{}.orig".format(f))
            os.chdir(self.parent_job_dir)                
            
    def run(self):
        """
        move to the job_dir, launch the job and back to the
         parent job directory
        """
        os.chdir(os.path.abspath(self.job_dir))
        logger.info('running in : '+self.job_dir)
        p = None
        #if launching jobs via batch system
        if self.vis.qadapter is not None:
            submit_cmd = self.vis.qadapter.supported_q_types[self.vis.qadapter.q_type]
            #print os.path.exists(self.vis.script_name)
            cmd = [submit_cmd, self.vis.script_name]
            with open(self.output_file, 'w') as f:            
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                                     stderr=subprocess.PIPE)
                stdout, stderr= p.communicate() 
                self.job_id = stdout.rstrip('\n')
                f.write(self.job_id)
            #reservation_id = self.vis.qadapter.\
            #submit_to_queue(self.vis.script_name)
            #cmd = ['echo', str(reservation_id)]
            #with open(self.output_file, 'w') as f:
            #    p = subprocess.Popen(cmd, stdout=f)
        else:
            cmd = list(self.job_cmd)
            with open(self.output_file, 'w') as f:
                p = subprocess.Popen(cmd, stdout=f, stderr=f)
            self.job_id = 0 #None
        os.chdir(self.parent_job_dir)
        if self.wait:
            return p
        else:
            return 0
    
    def postprocess(self):
        pass

    def name(self):
         return self.__class__.__name__

class MPINTVaspErrors(ErrorHandler):
    """
    handles restarting of jobs that exceed the walltime
    employs the check + correct method of custodian ErrorHandler
    """
    def check():
        """
        method to check completion of the job
        """
        pass

    def correct():
        """
        method to correct the error or restart the job
        """
        pass 
