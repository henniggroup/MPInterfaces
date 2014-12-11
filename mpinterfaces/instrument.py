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
from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vaspio_set import DictVaspInputSet #MPGGAVaspInputSet
from custodian.custodian import Job, gzip_dir
from custodian.vasp.interpreter import VaspModder


VASP_INPUT_FILES = {"INCAR", "POSCAR", "POTCAR", "KPOINTS"}

VASP_OUTPUT_FILES = ['DOSCAR', 'INCAR', 'KPOINTS', 'POSCAR', 'PROCAR',
                     'vasprun.xml', 'CHGCAR', 'CHG', 'EIGENVAL', 'OSZICAR',
                     'WAVECAR', 'CONTCAR', 'IBZKPT', 'OUTCAR']


class MPINTVaspInputSet(DictVaspInputSet):
    """
    defines the set of input required for a vasp job i.e
    create INCAR, POSCAR, POTCAR & KPOINTS files
    
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
        d = super(MPINTVaspInputSet, self).as_dict()
        return d



class MPINTVaspJob(Job):
    """
    defines a vasp job i.e setup the required input files and lanuch the job
    
    Args:
       job_cmd : a list, the command to be issued in each job_dir eg: ['qsub', 'submit_job']
       setup_dir : directory that has the setup files for creating the rest of the vasp inputs
       job_dir : the directory from which the jobs will be launched
    
    """
    def __init__(self, job_cmd, output_file="job.out", setup_dir='.', job_dir='untitled', suffix="",
                 final=True, gzipped=False, backup=True,
                 vis=None, auto_npar=True,
                 auto_gamma=True, settings_override=None,
                 gamma_vasp_cmd=None, copy_magmom=False):

        self.job_cmd = job_cmd
        self.launch_cmd = self.job_cmd[0]
        if len(self.job_cmd) > 1:
            self.launch_script = self.job_cmd[1]        
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
        #if launching jobs via batch system
        if self.launch_cmd == 'qsub':
            script_file = '../'+self.launch_script
            if os.path.isfile(script_file):
                shutil.copy(script_file, '.')
            else:
                print 'the submit script doesnt exit'
                sys.exit()
        cmd = list(self.job_cmd)
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

     
#test
#if __name__ == '__main__':

