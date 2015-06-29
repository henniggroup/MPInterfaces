from __future__ import division, unicode_literals, print_function

"""
Calibration module:  

This module contains the classes for 
1. Calibrate: Base class for specifying the parameters for 
calibration and setting up the VASP jobs in directory 
structure according to 
2. CalibrateBulk: calibrating a periodic bulk structure,
3. CalibrateSlab: creates a slab of given crystallographic facet, 
thickness and vacuum spacing, 
3. CalibrateMolecule: creates a molecule in a box 
4. CalibrateInterface: calibrates an interface composed of slab plus 
molecule

The attribute turn_knobs controls the parameters to be calibrated
for a given structure 

"""

import sys
import os
import shutil
import subprocess
import operator
import re
import time
import datetime
from itertools import product
from collections import Counter, OrderedDict
import logging

import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.io.vaspio.vasp_input import Incar, Poscar
from pymatgen.io.vaspio.vasp_input import Potcar, Kpoints
from pymatgen.io.vaspio.vasp_output import Outcar
from pymatgen.apps.borg.queen import BorgQueen

from custodian.vasp.handlers import VaspErrorHandler
#from custodian.vasp.handlers import FrozenJobErrorHandler
#from custodian.vasp.handlers import MeshSymmetryErrorHandler
#from custodian.vasp.handlers import NonConvergingErrorHandler
from custodian.custodian import Custodian, gzip_dir
from custodian.vasp.interpreter import VaspModder

from mpinterfaces.instrument import MPINTVaspInputSet, MPINTVaspJob
from mpinterfaces.data_processor import MPINTVaspDrone
from mpinterfaces.interface import Interface, Ligand
from mpinterfaces.utils import get_ase_slab

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


class Calibrate(object):    
    """
    The base class for creating vasp work flows for
    calibrating the input parameters for different systems
    
    """
    def __init__(self, incar, poscar, potcar, kpoints, system=None,
                 is_matrix = False, Grid_type = 'A',
                 setup_dir='.', parent_job_dir='.',job_dir='Job',
                 qadapter=None, job_cmd='qsub', wait=True,
                 turn_knobs=OrderedDict( [ ('ENCUT',[]),
                                           ('KPOINTS',[])] ) ):
        """
        Calibrate constructor

        Args:
            incar (Incar object): input INCAR
            poscar (Poscar object): input POSCAR
            potcar (Potcar object): input POTCAR
            kpoints: input KPOINTS
            system: system info as a dictionary,
                slab or interface example:
                system={'hkl':[1,1,1], 'ligand':None}, 
            is_matrix (bool): whether the jobs are dependent on each
                other
            Grid_type: kpoints grid_type
            setup_dir: directory from which the intial setup files
                will be copied from
            parent_job_dir: the directory from which all the jobs are
                launched
            job_dir: job directory
            qadapter: adapter for the batch system
            job_cmd: command to be used for submitting the job. If 
                qadapter is specified then job_cmd is ignored
            wait: whther to wait for the job to finish. If the job is
                being submitted to the queue then there is no need for
                waiting
            turn_knobs: an ordered dictionary of parmaters and the 
                corresponding values

        Note: input structure if needed will be obtained from the
            provided poscar object
        """
        self.name = datetime.datetime.now().isoformat()
        self.system = system
        self.parent_job_dir = os.path.abspath(parent_job_dir)
        self.setup_dir = os.path.abspath(setup_dir)
        self.job_dir = job_dir
        self.incar = incar               
        self.poscar = poscar
        self.potcar = Potcar(poscar.site_symbols)
        self.kpoints = kpoints
        self.incar_orig = incar.as_dict()
        self.poscar_orig = poscar.as_dict()
        self.potcar_orig = potcar.as_dict()
        self.kpoints_orig = kpoints.as_dict()
        self.qadapter = qadapter
        self.job_dir_list = []
        self.jobs = []
        self.job_ids = []        
        #example:- handlers = [VaspErrorHandler(),
        #FrozenJobErrorHandler(),
        #MeshSymmetryErrorHandler(), NonConvergingErrorHandler()]
        self.handlers = [] #[VaspErrorHandler("stdout")]
        self.job_cmd = job_cmd
        self.n_atoms = 0
        self.turn_knobs = turn_knobs
        self.response_to_knobs = {}
        self.sorted_response_to_knobs = {}        
        for k, v in turn_knobs.items():
            self.response_to_knobs[k] = {}
            self.sorted_response_to_knobs[k] = {}            
        self.is_matrix = is_matrix
        self.Grid_type = Grid_type
        self.wait = wait
    
    def setup(self):
        """
        set up the jobs for the given turn_knobs dict
        is_matrix = True implies that the params in the dict are 
        interrelated. Otherwise calcs corresponding to each dict key
        is independent
        """
        if self.is_matrix: 
            self.setup_matrix_job()	
        else:
            self._setup()

    def _setup(self, turn_knobs=None):
        """
        invoke the set up methods corresponding to the dict keys
        any key other than KPOINTS, VOLUME and POTCAR are treated
        as INCAR parameters
        
        Args:
            turn_knobs: knobs aka paramters to be tuned

        Note: poscar jobs setup through the VOLUME is only for
              backward compatibility, use POSCAR key in the
              turn_knobs to tune poscars
        """
        if turn_knobs is None:
            turn_knobs = self.turn_knobs
        if any(turn_knobs.values()):
            for k, v in turn_knobs.items():
                if k == 'KPOINTS' and v:
                    self.setup_kpoints_jobs(kpoints_list = v)
                elif k == 'VOLUME' and v:
                    self.setup_poscar_jobs(scale_list = v)
                elif k == 'POTCAR' and v:
                    self.setup_potcar_jobs(mappings = v)
                elif k == 'POSCAR' and v:
                    self.setup_poscar_jobs(poscar_list=v)
                else:
                    self.setup_incar_jobs(k, v)
        else:
            logger.warn('knobs not set, running a single job')
            self.add_job(name='single_job', job_dir=self.job_dir)

    def setup_matrix_job(self):
        """
        set up jobs where the dict keys are interrelated
        mind: its an ordered dict, the order in which the keys
        are specified determines the nested directory structure
        """
        orig_job_dir = self.job_dir
        job_dir = self.job_dir
        n_items = len(self.turn_knobs.items())
        keys = self.turn_knobs.keys()
        self._setup(turn_knobs=dict([(keys[0], 
                                      self.turn_knobs[keys[0]])]))
        self.recursive_jobs(n_items, keys, 0)
        #restore
        self.job_dir = orig_job_dir

    def recursive_jobs(self,n, keys, i):
        """
        recursively setup the jobs: used by setup_matrix_job
        
        Args:
            n: total number of knobs aka parameters to be tuned
            keys: list of knobs i.e parameter names
            i: ith knob
            
        """
        job_dir = self.job_dir + os.sep + self.key_to_name(keys[i])
        if i == n-1 and i != 0:
            for val in self.turn_knobs[keys[i]]:
                self.job_dir = job_dir + os.sep + self.val_to_name(val)
                logger.info('setting jobs in the directory: '+self.job_dir)
                self._setup(turn_knobs=dict([(keys[i], [val])]))            
                self.add_job(name=job_dir, job_dir=self.job_dir)
        else:
            for val in self.turn_knobs[keys[i]]:
                self.job_dir = job_dir + os.sep + self.val_to_name(val)
                logger.info('setting jobs in the directory: '+self.job_dir)
                self._setup(turn_knobs=dict([(keys[i], [val])]))
                self.recursive_jobs(n,keys,i+1)

    def key_to_name(self, key):
        """
        convenient string mapping for the keys in the turn_knobs dict
        
        Args:
            key: key to the knob dict
            
        Returns:
            an appropriate string representation of the key so that
            the name doesnt clash with the filenames
        """
        if key == 'KPOINTS':
            return 'KPTS'
        elif key == 'POTCAR':
            return 'POT'
        elif key == 'POSCAR':
            return 'POS'        
        else:
            return key

    def val_to_name(self, val):
        """
        convert a value to a string so that it can be used for naming
        the job directory
        the decimal points in floats are replaced with underscore 
        character
        if the value is of type list, kpoint_to_name method is used 
        since
        only kpoint values are expected to be of type list
        if the values is of type dict then potcar_to_name method is 
        invoked
        
        Args:
            val: knob value to be converted into an appropriate string
                representation
            
        Returns:
            a string filename for the value 
        """
        if type(val) == float:
            return re.sub('\.','_',str(val))
        elif type(val) == list:
            return self.kpoint_to_name(val, 'M')
        elif type(val) == dict:
            return self.potcar_to_name(val)
        elif isinstance(val, Poscar):
            return val.comment        
        else:
            return str(val)
                
    def kpoint_to_name(self, kpoint, grid_type):
        """
        get a string representation for the given kpoint
        
        Args:
            kpoint: an iterable
            grid_type: grid_type used for the KPOINTS
        
        Returns:    
            string representation for kpoint eg: Monkhorst Pack
	    2 2 2 will be named 2x2x2
        """
        if grid_type == 'M':
            return str(kpoint[0]) + 'x' + str(kpoint[1]) + 'x' \
                + str(kpoint[2])
        elif grid_type == 'A':    
            return str(kpoint)

    def potcar_to_name(self, mapping):
        """
        convert a symbol mapping to a name that can be used for
        setting up the potcar jobs
         
         Args:
             mapping: example:- if mapping = {'Pt':'Pt_pv', 
                 'Si':'Si_GW'} then the name will be Pt_pvSi_GW
                 
        Returns:
            string 
        """
        l = [v for k,v in mapping.items()]
        return ''.join(l)
        

    def set_incar(self, param, val):
        """
        set the incar paramter, param = val
        """
        self.incar[param] = val

    def set_poscar(self, scale=None, poscar=None):
        """
        perturbs given structure by volume scaling factor 
	or takes user defined variants of Poscar 
	
	Args:
	   scale : Volume Scaling parameter

	   poscar : Poscar object of user defined structure

	   set the poscar: volume scaled by the scale factor
        """
        if scale is not None:
            structure = Poscar.from_dict(self.poscar_orig).structure
            volume = structure.volume
            structure.scale_lattice(scale *volume)
            self.poscar = Poscar(structure)
        elif poscar is not None:
            self.poscar = poscar

    def set_potcar(self, mapping):
        """
        set the potcar: symbol to potcar type mapping
        """
        symbols = self.poscar.site_symbols
        mapped_symbols = []
        for sym in symbols:
            mapped_symbols.append(mapping[sym])
        self.potcar = Potcar(symbols=mapped_symbols)
        pass

    def set_kpoints(self, kpoint):
        """
        set the kpoint
        """
        if self.Grid_type == 'M':
            self.kpoints = Kpoints.monkhorst_automatic(kpts = kpoint)
        elif self.Grid_type == 'A':
            self.kpoints = Kpoints.automatic(subdivisions = kpoint)
        name = self.kpoint_to_name(kpoint, self.Grid_type)
        logger.info('KPOINTS = '+name)
        job_dir = self.job_dir +os.sep+ self.key_to_name('KPOINTS') \
          + os.sep + name
        return job_dir
                                        
    def setup_incar_jobs(self, param, val_list):
        """
        set up incar jobs, called through a dictionary of {param: val_list} 
	by turn_knobs, calls set_incar to set the value to param

	Args:
	    param: Name of INCAR parameter
	    val_list: List of values to vary for the param
        """
        if val_list:
            for val in val_list:
                logger.info('setting INCAR parameter ' + param + ' = '\
                            + str(val))
                self.set_incar(param, val)
                if not self.is_matrix:
                    job_dir  = self.job_dir+ os.sep + \
                        param + os.sep +  self.val_to_name(val)
                    self.add_job(name=job_dir, job_dir=job_dir)
        else:
        	logger.warn('incar list empty')
                    
            
    def setup_kpoints_jobs(self, kpoints_list = []):
        """
        setup the kpoint jobs
        
        """
        if kpoints_list:
            for kpoint in kpoints_list:
                job_dir = self.set_kpoints(kpoint)
                if not self.is_matrix:                     
                    self.add_job(name=job_dir, job_dir=job_dir)
        else:
        	logger.warn('kpoints_list empty')
            
    def setup_poscar_jobs(self, scale_list=[], poscar_list=[]):
        """
        for scaling the latice vectors of the original structure,
        scale_list is volume scaling factor list
        """
        if scale_list:
            for scale in scale_list:
                self.set_poscar(scale=scale)
                job_dir  = self.job_dir+ os.sep + 'POS' +\
                        os.sep + 'VOLUME_'+str(scale)
                if not self.is_matrix:
                    self.add_job(name=job_dir, job_dir=job_dir)
        elif poscar_list:
            for poscar in poscar_list:
                self.set_poscar(poscar=poscar)
                job_dir  = self.job_dir+ os.sep +'POS' +\
                  os.sep + poscar.comment
                if not self.is_matrix:
                    self.add_job(name=job_dir, job_dir=job_dir)
                    

    def setup_potcar_jobs(self, mappings):
        """
        take a list of symbol mappings and setup the potcar jobs
        """
        for mapping in mappings:
                self.set_potcar(mapping)
                if not self.is_matrix:
                    job_dir  = self.job_dir+ os.sep \
                        + self.key_to_name('POTCAR') \
                        + os.sep + self.potcar_to_name(mapping)
                    self.add_job(name=job_dir, job_dir=job_dir)
                

    def add_job(self, name='noname', job_dir='.'):
        """
        add a single job using the current incar, poscar, potcar and 
        kpoints
        """
        vis = MPINTVaspInputSet(name, self.incar, self.poscar,
                                self.potcar, self.kpoints,
                                self.qadapter)
        #the job command can be overrridden in the run method
        job = MPINTVaspJob(self.job_cmd, name=name, final = True,
                           setup_dir=self.setup_dir,
                           parent_job_dir=self.parent_job_dir,
                           job_dir=job_dir, vis=vis, auto_npar=False,
                            auto_gamma=False, wait=self.wait)
        self.job_dir_list.append(os.path.abspath(job_dir))
        self.jobs.append(job)
    
    def run(self, job_cmd=None):
        """
        run the vasp jobs through custodian
        if the job list is empty,
        run a single job with the initial input set
        """
        for j in self.jobs:
            if job_cmd is not None:            
                j.job_cmd = job_cmd
            else:
                j.job_cmd = self.job_cmd
        c_params = {'jobs': [j.as_dict() for j in self.jobs],
                'handlers': [h.as_dict() for h in self.handlers],
                'max_errors': 5}
        c = Custodian(self.handlers, self.jobs, max_errors=5)
        c.run()
        for j in self.jobs:
            self.job_ids.append(j.job_id)

    def set_knob_responses(self):
        """
        set up a dictionary that maps the turn knob keys and
        their values to the calculated energy
        """
        drone = \
          MPINTVaspDrone(inc_structure=True, inc_incar_n_kpoints=True)
        bg =  BorgQueen(drone)
        for k, v in self.response_to_knobs.items():
            rootpath = self.job_dir+ os.sep + k
            logger.info('rootpath = '+rootpath)
            #bg.parallel_assimilate(rootpath)        
            bg.serial_assimilate(rootpath)
            allentries =  bg.get_data()
            for e in allentries:
                if e:
                    self.n_atoms = len(e.structure)
                    if k == 'KPOINTS':
                        self.response_to_knobs[k][str(e.kpoints.kpts)] \
                           = e.energy
                    else:
                        self.response_to_knobs[k][str(e.incar[k])] \
                          = e.energy

    def set_sorted_optimum_params(self):
        """
        sort the dictionary of energy values and enforce the
        convergence criterion.
        Finally, get the optimum parameter values from the set of
        values that satisfy the convergence criterion
        """
        matching_knob_responses = []
        sorted_knob_responses = []                
        self.optimum_knob_responses = {}
        #order the keys(encut or kpoint)from large value of the
        # energy  to small value
        for k, v in self.response_to_knobs.items():
            sorted_knob_responses = \
              sorted(v.items(), 
                     key=operator.itemgetter(1), reverse=True)
            #get the list of encut and kpoints that
            #satisfy the delate criterion
            #mind: default deltae = 0.001eV per atom
            matching_knob_responses \
                = self.enforce_cutoff(sorted_knob_responses)
            self.sorted_response_to_knobs[k] \
                = OrderedDict(sorted_knob_responses)
            if matching_knob_responses:
                if k == "KPOINTS" and self.Grid_type == 'M':
                    nkpt = matching_knob_responses[0][0] * \
                      matching_knob_responses[0][1] * matching_knob_responses[0][2]
                    for i, val in enumerate(matching_knob_responses):
                        if i < len(matching_kpt)-1:
                            nkpt1 = val[i+1][0] * val[i+1][1] * val[i+1][2]
                            if nkpt1<nktp:
                                self.optimum_knob_responses[k] = val
                else:
                    self.optimum_knob_responses[k] = min(matching_knob_responses)
                        
    def enforce_cutoff(self, input_list, delta_e_peratom=0.001):
        """
        enfore convergence criterion: energy difference of 1meV per
        atom.
        
        returns a list of the parameters that satisfy the criterion
        """
        matching_list = []
        for i, e in enumerate(input_list):
            if i < len(input_list)-1:
                if np.abs(input_list[i+1][1] - e[1])/self.n_atoms <= \
                  delta_e_peratom:
                    matching_list.append(input_list[i+1][0])
        if matching_list:
            matching_kpt_list = []
            if '[[' in matching_list[0]:
                for ml in matching_list:
                    if '[[' in ml:
                        m = re.search(r"\[\[(\d+)\,(\d+)\,(\d+)\]\]",
                                       ml)
                        matching_kpt_list.append( [ int(m.group(1)),
                                                    int(m.group(2)),
                                                     int(m.group(3))])
                return matching_kpt_list
            else:
                return [float(val) for val in matching_list]
        else:
            return []
        
    @staticmethod
    def check_calcs(cal_objs):
        """
        checks the OUTCAR file to see whether the calulation is done
        or not. Also checks(rather naively) for whether the
        calculation is running or not.
        Sets the calc_done and isrunning variables in the calibrate 
        objects
        """
        done = []
        i = 0
        for cal in cal_objs:        
            cal.calc_done = False
            cal.isrunning = False
            for jdir in cal.job_dir_list:
                i += 1
                outcar_file = jdir + os.sep + 'OUTCAR'
                if os.path.isfile(outcar_file):
                    mtime = os.stat(outcar_file)[8]
                    last_mod_time =  datetime.datetime.fromtimestamp(mtime)
                    current_time = datetime.datetime.now()
                    logger.info('time delta {}'
                                .format(current_time - last_mod_time))
                    #check whether the OUTCAR file had been modified in the
                    # last hour
                    #if it had not been modified in the past hour, the 
                    #calculation is assumed dead
                    if current_time - last_mod_time < datetime.timedelta(seconds=3600):
                        cal.isrunning = True
                    outcar = Outcar(outcar_file)
                    for k in outcar.run_stats.keys():
                        if 'time' in k:
                            cal.calc_done = True
                            done.append(cal.calc_done)
                            logger.info('done in {}'.format(jdir))
                            break
        if done:
            if len(done) == i:
                return all(done)
            else:
                return False
        else:
            return False
    
    def as_dict(self):
        d = {}
        d['calibrate'] = self.__class__.__name__
        d['name'] = self.name
        d['incar'] = self.incar_orig
        d['poscar'] = self.poscar_orig
        d['kpoints'] = self.kpoints_orig
        d['turn_knobs'] = self.turn_knobs
        d['job_dir_list'] = self.job_dir_list
        d['job_ids'] = self.job_ids        
        if self.system is not None:
            d['system'] = self.system
        return d
        
class CalibrateMolecule(Calibrate):
    """
    
    Calibrate paramters for Molecule calculations
    
    """
    def __init__(self, incar, poscar, potcar, kpoints, system=None,
                 is_matrix = False, Grid_type = 'A',
                 setup_dir='.', parent_job_dir='.',
                 job_dir='./Molecule', qadapter=None,
                 job_cmd='qsub', wait=True,
                 turn_knobs={'ENCUT':[],'KPOINTS':[]}):
        
        Calibrate.__init__(self, incar, poscar, potcar, kpoints, 
                           system=system, is_matrix = is_matrix, 
                           Grid_type = Grid_type, setup_dir=setup_dir,
                           parent_job_dir=parent_job_dir,
                           job_dir=job_dir, qadapter=qadapter,
                           job_cmd=job_cmd, wait=wait,
                           turn_knobs = turn_knobs)
        
    def setup_kpoints_jobs(self, Grid_type = 'M',
                           kpoints_list = None, conv_step = 1):
        logger.warn("Its a molecule ! no need for kpoint convergence")
        self.kpoints = Kpoints.monkhorst_automatic(kpts = [1,1,1])
        return

    
class CalibrateBulk(Calibrate):
    """
    
    Calibrate parameters for Bulk calculations
    
    """
    def __init__(self, incar, poscar, potcar, kpoints, system=None,
                 is_matrix = False, Grid_type = 'A',
                  setup_dir='.', parent_job_dir='.',
                  job_dir='./Bulk', qadapter=None,
                  job_cmd='qsub', wait=True,
                  turn_knobs={'ENCUT':[],'KPOINTS':[]}): 
            
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                           system=system, is_matrix = is_matrix,
                           Grid_type = Grid_type, setup_dir=setup_dir,
                           parent_job_dir=parent_job_dir,
                           job_dir=job_dir, qadapter=qadapter,
                           job_cmd=job_cmd,wait=wait,
                           turn_knobs = OrderedDict(turn_knobs))
        

class CalibrateSlab(Calibrate):
    """
    Calibrate paramters for Slab calculations
    """
    def __init__(self, incar, poscar, potcar, kpoints, system=None,
                 is_matrix = False, Grid_type = 'A',
                 setup_dir='.', parent_job_dir='.', job_dir='./Slab',
                 qadapter=None, job_cmd='qsub', wait=True,
                 turn_knobs={'VACUUM':[],'THICKNESS':[]}, from_ase=False):
        self.from_ase = from_ase
        self.is_matrix = is_matrix
        self.system = system
        self.input_structure = poscar.structure.copy()
        self.slab_setup(turn_knobs=turn_knobs)
        Calibrate.__init__(self, incar, poscar, potcar, kpoints, 
                           system=system, is_matrix = is_matrix,
                           Grid_type = Grid_type,setup_dir=setup_dir,
                           parent_job_dir=parent_job_dir,
                           job_dir=job_dir, qadapter=qadapter,
                           job_cmd=job_cmd, wait=wait,
                           turn_knobs = turn_knobs)

    def slab_setup(self, turn_knobs=None):
        """
        invoke the set up methods corresponding to the dict keys:
        VACUUM and THICKNESS
        sets the POSCAR key in the turn_knobs 
        """
        if turn_knobs is None:
            turn_knobs = self.turn_knobs
        if any(turn_knobs.values()):
            keys = ['VACUUM', 'THICKNESS']
            poscar_list = []
            if self.is_matrix:
                prod_list = [turn_knobs[k] for k in keys]
                for params in product(*tuple(prod_list)):
                    poscar_list.append(self.create_slab(*params))
            else:                    
                for k, v in turn_knobs.items():
                    if k == 'VACUUM' and v:
                        poscar_list += self.setup_vacuum_jobs(v)
                    elif k == 'THICKNESS' and v:
                        poscar_list += self.setup_thickness_jobs(v)
            for k in keys:
                if turn_knobs.get(k):
                    del turn_knobs[k]
            turn_knobs['POSCAR'] = poscar_list

    def setup_vacuum_jobs(self, vacuum_list):
        """
        create slabs with the provided vacuum settings

        returns list of poscars        
        """
        return [self.create_slab(vacuum=val) for val in vacuum_list]

    def setup_thickness_jobs(self, thickness_list):
        """
        create slabs with the provided thickness settings

        returns list of poscars
        """
        return [self.create_slab(thickness=val) for val in thickness_list]
            
    def create_slab(self, vacuum=12, thickness=10):
        """
        set the vacuum spacing, slab thickness and call sd_flags
        for top 2 layers

        returns the poscar corresponding to the modified structure
        """ 
        strt_structure = self.input_structure.copy()
        if self.from_ase:
            slab_struct = get_ase_slab(strt_structure, hkl=self.system['hkl'],
                                        min_thick=thickness, min_vac=vacuum)
        else:
            slab_struct= SlabGenerator(initial_structure= strt_structure,
                                       miller_index= self.system['hkl'], 
                                       min_slab_size= thickness,
                                       min_vacuum_size=vacuum, 
                                       lll_reduce=False, center_slab=True,
                                       primitive=False).get_slab()
        slab_struct.sort()
        sd = self.set_sd_flags(slab_struct)
        comment = 'VAC'+str(vacuum)+'THICK'+str(thickness)
        return Poscar(slab_struct, comment=comment,
                      selective_dynamics=sd)    

    @staticmethod
    def set_sd_flags(interface=None, n_layers=2):
        """
        set the relaxation flags for top and bottom layers of interface.
        
        The upper and lower bounds of the z coordinate are determined
        based on the slab. All layers above and below the bounds will
        be relaxed. This means if there is a ligand on top of the slab,
        all of its atoms will also be relaxed.
        """
        sd_flags = np.zeros_like(interface.frac_coords)
        if isinstance(interface, Interface):
            slab = interface.slab
        else:
            slab = interface
        z_coords = slab.frac_coords[:,2]
        z_lower_bound = np.unique(z_coords)[n_layers-1]
        z_upper_bound = np.unique(z_coords)[-n_layers]
        sd_flags[ [i for i, coords in enumerate(slab.frac_coords)
                  if coords[2]>=z_upper_bound or coords[2]<=z_lower_bound] ] \
                  = np.ones((1,3))
        return sd_flags

    def set_reconstructed_surface(self, sites_to_add):
        """
        Append sites as needed for reconstruction TODO

        """
        pass  

    
class CalibrateInterface(CalibrateSlab):
    """
    
    Calibrate paramters for interface calculations
    
    """
    def __init__(self, incar, poscar, potcar, kpoints, system=None,
                 is_matrix = False, Grid_type = 'A',
                 setup_dir='.', parent_job_dir='.', job_dir='./Interface',
                 qadapter=None, job_cmd='qsub', wait=True,
                 turn_knobs={'VACUUM':[],'THICKNESS':[]}, from_ase=False):
        CalibrateSlab.__init__(self, incar, poscar, potcar, kpoints, 
                           system=system, is_matrix = is_matrix, 
                           Grid_type = Grid_type, setup_dir=setup_dir,
                           parent_job_dir=parent_job_dir,
                           job_dir=job_dir, qadapter=qadapter,
                           job_cmd=job_cmd, wait=wait,
                           turn_knobs = turn_knobs, from_ase=from_ase)
        self.interface_setup(turn_knobs=turn_knobs)        

    def interface_setup(self, turn_knobs=None):
        if self.system['ligand'] is None:
            return
        else:
            if turn_knobs is None:
                turn_knobs = self.turn_knobs
            if any(turn_knobs.values()):
                poscar_list = []
                poscar_list.append(self.create_interface())
                turn_knobs['POSCAR'] = poscar_list
            
    def create_interface(self):
        """
        add params that you want to vary
        """
        structure = self.input_structure.copy()
        iface = Interface(structure,
                          hkl=self.system['hkl'],
                          ligand = Ligand.from_dict(self.system['ligand']),
                          from_ase=self.from_ase)
        iface.sort()
        sd = self.set_sd_flags(iface, n_layers=2)
        #if theer are other paramters that are being varied
        #change the comment accordingly
        comment = self.system['hkl']+self.system['ligand']['name']
        return Poscar(slab_struct, comment=comment,
                      selective_dynamics=sd)
            
if __name__ == '__main__':
    #STRUCTURE
    a0 = 3.965
    lvec = [ [0.5, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5] ]
    lvec = np.array(lvec) * a0
    lattice = Lattice(lvec)
    structure = Structure( lattice, ['Pt'], [ [0.0, 0.0, 0.0] ],
                           coords_are_cartesian=False )

    #INITIAL VASP INPUT SET
    incarparams = {'System':'test',
                   'ENCUT': 400,
                   'ISMEAR': 1,
                   'SIGMA': 0.1,
                   'EDIFF': 1E-6}
    incar = Incar(params = incarparams)
    poscar = Poscar(structure, comment='test')
    potcar = Potcar( symbols = poscar.site_symbols, functional='PBE',
                     sym_potcar_map=None)
    kpoints = Kpoints.monkhorst_automatic(kpts=(16, 16, 16),
                                          shift=(0, 0, 0))

    #CALIBRATION INPUT
    system = {'hkl':[1,1,1], 'ligand':None }    
    turn_knobs = OrderedDict( [
        ('SIGMA', [0.025, 0.50]),
        ('POTCAR', [{'Pt':'Pt'}, {'Pt':'Pt_pv'}, {'Pt':'Pt_GW'}]),
        ('IBRION', [1, 2, 3]),
        ('KPOINTS', [k for k in range(20, 40, 10)]),
        ('ENCUT', range(400,700,100)),
        ('VACUUM', [10,12,15]),
        ('THICKNESS', [11])        
        ] )
    is_matrix = True
    job_dir = 'Slab'
    job_cmd = ['ls','-lt'] 

    #SETUP AND RUN    
    cal = CalibrateSlab( incar, poscar, potcar, kpoints,
                         system = system,
                         is_matrix = is_matrix,
                         job_dir = job_dir,
                         turn_knobs = turn_knobs )
    cal.setup()     
    cal.run(job_cmd)
    #print(cal.job_ids)

    
    # to use the next after the calibrate jobs are done
    #done = Calibrate.check_calcs([cal])  #
    #get the knob responses
    #cal.set_knob_responses()
    #optimu knob responses
    #cal.set_sorted_optimum_params()
    #print(cal.sorted_response_to_knobs['ENCUT']['600.0'])
    #print(cal.optimum_knob_responses)
    #test enforce_cutoff
    #setup relaxation jobs after optimum parameters are set
    #cal.setup_relaxation_jobs([cal])
    #inp_list = [ ['[[2,2,4]]', 10], ['[[2,2,5]]', 9.9],
    #           ['[[2,2,6]]', 9.895], ['[[2,2,7]]', 9.888],
    #['[[2,2,8]]', 9.879],]
    #print(cal.enforce_cutoff(inp_list, delta_e=0.01))
