# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

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
import re
import datetime
from itertools import product
from collections import OrderedDict
import logging

import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.inputs import Potcar, Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath

from custodian.custodian import Custodian

from monty.json import MSONable, MontyEncoder
from monty.serialization import dumpfn

from mpinterfaces.instrument import MPINTVaspInputSet, MPINTVaspJob
from mpinterfaces.interface import Interface, Ligand
from mpinterfaces.utils import get_ase_slab

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


class Calibrate(MSONable):
    """
    The base class for creating vasp work flows for
    calibrating the input parameters for different systems

    A wrapper around Custodian
    """
    LOG_FILE = "calibrate.json"

    def __init__(self, incar, poscar, potcar, kpoints, system=None,
                 is_matrix=False, Grid_type='A',
                 parent_job_dir='.', job_dir='Job',
                 qadapter=None, job_cmd='qsub', wait=True,
                 mappings_override=None, functional="PBE",
                 turn_knobs=OrderedDict([('ENCUT', []),
                                         ('KPOINTS', [])]),
                 checkpoint_file=None, cal_logger=None):
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
            mappings_override: override symbol mapping in potcar
                               eg:- {'S':'S_sv'}
            functional: exchange-correlation functional

        Note: input structure if needed will be obtained from the
            provided poscar object
        """
        self.name = datetime.datetime.now().isoformat()
        self.system = system
        self.parent_job_dir = os.path.abspath(parent_job_dir)
        self.job_dir = job_dir
        self.incar = incar
        self.poscar = poscar
        self.potcar = potcar
        if poscar:
            self.potcar = Potcar(symbols=poscar.site_symbols,
                                 functional=functional)
        self.kpoints = kpoints
        if incar:
            self.incar_orig = incar.as_dict()
        if poscar:
            self.poscar_orig = poscar.as_dict()
        if self.potcar:
            self.potcar_orig = self.potcar.as_dict()
        if kpoints:
            self.kpoints_orig = kpoints.as_dict()
        self.qadapter = qadapter
        self.job_dir_list = []
        self.jobs = []
        self.job_ids = []
        self.handlers = []
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
        self.cal_log = []
        self.mappings_override = mappings_override
        self.functional = functional
        self.checkpoint_file = checkpoint_file
        if cal_logger:
            self.logger = cal_logger
        else:
            self.logger = logger

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
                    self.setup_kpoints_jobs(kpoints_list=v)
                elif k == 'VOLUME' and v:
                    self.setup_poscar_jobs(scale_list=v)
                elif k == 'POTCAR' and v:
                    self.setup_potcar_jobs(mappings=v)
                elif k == 'POSCAR' and v:
                    self.setup_poscar_jobs(poscar_list=v)
                else:
                    self.setup_incar_jobs(k, v)
        else:
            self.logger.warn('knobs not set, running a single job')
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
        # restore
        self.job_dir = orig_job_dir

    def recursive_jobs(self, n, keys, i):
        """
        recursively setup the jobs: used by setup_matrix_job
        
        Args:
            n: total number of knobs aka parameters to be tuned
            keys: list of knobs i.e parameter names
            i: ith knob
            
        """
        job_dir = self.job_dir + os.sep + self.key_to_name(keys[i])
        if i == n - 1 and i != 0:
            for val in self.turn_knobs[keys[i]]:
                self.job_dir = job_dir + os.sep + self.val_to_name(val)
                self.logger.info('setting jobs in the directory: ' + self.job_dir)
                self._setup(turn_knobs=dict([(keys[i], [val])]))
                self.add_job(name=job_dir, job_dir=self.job_dir)
        else:
            for val in self.turn_knobs[keys[i]]:
                self.job_dir = job_dir + os.sep + self.val_to_name(val)
                self.logger.info('setting jobs in the directory: ' + self.job_dir)
                self._setup(turn_knobs=dict([(keys[i], [val])]))
                self.recursive_jobs(n, keys, i + 1)

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
        elif key == 'POTCAR_map' or key == 'POTCAR_functional':
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
            return re.sub('\.', '_', str(val))
        elif type(val) == list:
            return self.kpoint_to_name(val, 'M')
        elif type(val) == dict:
            return self.potcar_to_name(val)
        elif isinstance(val, Poscar):
            return str(val.structure.composition.reduced_formula) \
                   + '_' + str(int(val.structure.lattice.volume)) \
                   + '_' + ''.join((val.comment).split())
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
        if grid_type == 'M' or grid_type == 'G':
            return str(kpoint[0]) + 'x' + str(kpoint[1]) + 'x' \
                   + str(kpoint[2])
        else:
            return str(kpoint)

    def potcar_to_name(self, mapping, functional):
        """
        convert a symbol mapping and functional to a name that 
        can be used for setting up the potcar jobs
         
         Args:
             mapping: example:- if mapping = {'Pt':'Pt_pv', 
                 'Si':'Si_GW'} then the name will be PBE_Pt_pv_Si_GW
                  with self.functional="PBE"
                 
        Returns:
            string 
        """
        if mapping:
            l = [v for k, v in mapping.items()]
            return '_'.join(self.functional, l)
        elif functional:
            return '_'.join(functional)
        else:
            return '_'.join(self.functional)

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
            structure.scale_lattice(scale * volume)
            self.poscar = Poscar(structure)
        elif poscar is not None:
            self.poscar = poscar

    def set_potcar(self, mapping=None, functional='PBE'):
        """
        set the potcar: symbol to potcar type mapping
        """
        symbols = self.poscar.site_symbols
        mapped_symbols = []
        if mapping:
            for sym in symbols:
                mapped_symbols.append(mapping[sym])
        elif self.mappings_override:
            for sym in symbols:
                if sym in self.mappings_override.keys():
                    mapped_symbols.append(self.mappings_override[sym])
                else:
                    mapped_symbols.append(sym)
        else:
            mapped_symbols = symbols
        if functional:
            func = functional
        else:
            func = self.functional
        self.potcar = Potcar(symbols=mapped_symbols,
                             functional=func)

    def set_kpoints(self, kpoint):
        """
        set the kpoint
        """
        if self.Grid_type == 'M':
            self.kpoints = Kpoints.monkhorst_automatic(kpts=kpoint)
        elif self.Grid_type == 'A':
            self.kpoints = Kpoints.automatic(subdivisions=kpoint)
        elif self.Grid_type == 'G':
            self.kpoints = Kpoints.gamma_automatic(kpts=kpoint)
        elif self.Grid_type == '3DD':
            self.kpoints = Kpoints.automatic_density_by_vol(structure= \
                                                                self.poscar.structure, kppvol=kpoint)
        elif self.Grid_type == 'band':
            self.kpoints = Kpoints.automatic_linemode(divisions=kpoint, \
                                                      ibz=HighSymmKpath(self.poscar.structure))

    def setup_incar_jobs(self, param, val_list):
        """
        set up incar jobs,, calls set_incar to set the value to param

        Args:
            param: Name of INCAR parameter
            val_list: List of values to vary for the param
        """
        if val_list:
            for val in val_list:
                self.logger.info('setting INCAR parameter ' + param + ' = ' \
                                 + str(val))
                self.set_incar(param, val)
                if not self.is_matrix:
                    job_dir = self.job_dir + os.sep + \
                              param + os.sep + self.val_to_name(val)
                    self.add_job(name=job_dir, job_dir=job_dir)
        else:
            self.logger.warn('incar list empty')

    def setup_kpoints_jobs(self, kpoints_list=None):
        """
        setup the kpoint jobs
        
        """
        if kpoints_list:
            for kpoint in kpoints_list:
                self.set_kpoints(kpoint)
                if not self.is_matrix:
                    job_dir = self.job_dir + os.sep + self.key_to_name('KPOINTS') \
                              + os.sep + self.kpoint_to_name(kpoint, self.Grid_type)
                    self.add_job(name=job_dir, job_dir=job_dir)
        else:
            self.logger.warn('kpoints_list empty')

    def setup_poscar_jobs(self, scale_list=None, poscar_list=None):
        """
        for scaling the latice vectors of the original structure,
        scale_list is volume scaling factor list
        """
        if scale_list:
            for scale in scale_list:
                self.set_poscar(scale=scale)
                self.set_potcar()
                if not self.is_matrix:
                    job_dir = self.job_dir + os.sep + 'POS' + \
                              os.sep + 'VOLUME_' + str(scale)
                    self.add_job(name=job_dir, job_dir=job_dir)
        elif poscar_list:
            for poscar in poscar_list:
                self.set_poscar(poscar=poscar)
                self.set_potcar()
                if not self.is_matrix:
                    job_dir = self.job_dir + os.sep + 'POS' + \
                              os.sep + self.val_to_name(poscar)
                    self.add_job(name=job_dir, job_dir=job_dir)

    def setup_potcar_jobs(self, mappings, functional_list):
        """
        take a list of symbol mappings and setup the potcar jobs
        """
        if functional_list:
            for func in functional_list:
                self.set_potcar(functional=func)
                if not self.is_matrix:
                    job_dir = self.job_dir + os.sep \
                              + self.key_to_name('POTCAR') \
                              + os.sep + self.potcar_to_name(func)
                    self.add_job(name=job_dir, job_dir=job_dir)

        elif mappings:
            for mapping in mappings:
                self.set_potcar(mapping)
                if not self.is_matrix:
                    job_dir = self.job_dir + os.sep \
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
                                self.qadapter, vis_logger=self.logger)
        # the job command can be overrridden in the run method
        job = MPINTVaspJob(self.job_cmd, name=name, final=True,
                           parent_job_dir=self.parent_job_dir,
                           job_dir=job_dir, vis=vis, wait=self.wait,
                           vjob_logger=self.logger)
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
            self.cal_log.append({"job": j.as_dict(),
                                 'job_id': j.job_id,
                                 "corrections": [],
                                 'final_energy': None})
            self.job_ids.append(j.job_id)
        if self.checkpoint_file:
            dumpfn(self.cal_log, self.checkpoint_file,
                   cls=MontyEncoder, indent=4)
        else:
            dumpfn(self.cal_log, Calibrate.LOG_FILE, cls=MontyEncoder,
                   indent=4)

    def as_dict(self):
        qadapter = None
        system = None
        if self.qadapter:
            qadapter = self.qadapter.to_dict()
        if self.system is not None:
            system = self.system
        d = dict(incar=self.incar.as_dict(),
                 poscar=self.poscar.as_dict(),
                 potcar=self.potcar.as_dict(),
                 kpoints=self.kpoints.as_dict(),
                 system=system, is_matrix=self.is_matrix,
                 Grid_type=self.Grid_type,
                 parent_job_dir=self.parent_job_dir,
                 job_dir=self.job_dir,
                 qadapter=qadapter, job_cmd=self.job_cmd,
                 wait=self.wait,
                 turn_knobs=self.turn_knobs,
                 job_dir_list=self.job_dir_list,
                 job_ids=self.job_ids)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        # d['calibrate'] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        incar = Incar.from_dict(d["incar"])
        poscar = Poscar.from_dict(d["poscar"])
        potcar = Potcar.from_dict(d["potcar"])
        kpoints = Kpoints.from_dict(d["kpoints"])
        cal = Calibrate(incar, poscar, potcar, kpoints,
                        system=d["system"], is_matrix=d["is_matrix"],
                        Grid_type=d["Grid_type"],
                        parent_job_dir=d["parent_job_dir"],
                        job_dir=d["job_dir"], qadapter=d.get("qadapter"),
                        job_cmd=d["job_cmd"], wait=d["wait"],
                        turn_knobs=d["turn_knobs"])
        cal.job_dir_list = d["job_dir_list"]
        cal.job_ids = d["job_ids"]
        return cal


class CalibrateMolecule(Calibrate):
    """
    
    Calibrate paramters for Molecule calculations
    
    """

    def __init__(self, incar, poscar, potcar, kpoints, system=None,
                 is_matrix=False, Grid_type='A',
                 parent_job_dir='.',
                 job_dir='./Molecule', qadapter=None,
                 job_cmd='qsub', wait=True,
                 mappings_override=None, functional="PBE",
                 turn_knobs={'ENCUT': [], 'KPOINTS': []},
                 checkpoint_file=None, cal_logger=None):
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                           system=system, is_matrix=is_matrix,
                           Grid_type=Grid_type,
                           parent_job_dir=parent_job_dir,
                           job_dir=job_dir, qadapter=qadapter,
                           job_cmd=job_cmd, wait=wait,
                           mappings_override=mappings_override,
                           functional=functional,
                           turn_knobs=turn_knobs,
                           checkpoint_file=checkpoint_file,
                           cal_logger=cal_logger)

    def setup_kpoints_jobs(self, Grid_type='M',
                           kpoints_list=None, conv_step=1):
        self.logger.warn("Its a molecule ! no need for kpoint convergence")
        self.kpoints = Kpoints.monkhorst_automatic(kpts=[1, 1, 1])
        return


class CalibrateBulk(Calibrate):
    """
    
    Calibrate parameters for Bulk calculations
    
    """

    def __init__(self, incar, poscar, potcar, kpoints, system=None,
                 is_matrix=False, Grid_type='A',
                 parent_job_dir='.',
                 job_dir='./Bulk', qadapter=None,
                 job_cmd='qsub', wait=True,
                 mappings_override=None, functional="PBE",
                 turn_knobs={'ENCUT': [], 'KPOINTS': []},
                 checkpoint_file=None, cal_logger=None):
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                           system=system, is_matrix=is_matrix,
                           Grid_type=Grid_type,
                           parent_job_dir=parent_job_dir,
                           job_dir=job_dir, qadapter=qadapter,
                           job_cmd=job_cmd, wait=wait,
                           mappings_override=mappings_override,
                           functional=functional,
                           turn_knobs=OrderedDict(turn_knobs),
                           checkpoint_file=checkpoint_file,
                           cal_logger=cal_logger)


class CalibrateSlab(Calibrate):
    """
    Calibrate paramters for Slab calculations
    """

    def __init__(self, incar, poscar, potcar, kpoints, system=None,
                 is_matrix=False, Grid_type='A',
                 parent_job_dir='.', job_dir='./Slab',
                 qadapter=None, job_cmd='qsub', wait=True,
                 mappings_override=None, functional="PBE",
                 turn_knobs={'VACUUM': [], 'THICKNESS': []},
                 from_ase=False, checkpoint_file=None,
                 cal_logger=None):
        self.from_ase = from_ase
        self.is_matrix = is_matrix
        self.system = system
        self.input_structure = poscar.structure.copy()
        self.slab_setup(turn_knobs=turn_knobs)
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                           system=system, is_matrix=is_matrix,
                           Grid_type=Grid_type,
                           parent_job_dir=parent_job_dir,
                           job_dir=job_dir, qadapter=qadapter,
                           job_cmd=job_cmd, wait=wait,
                           mappings_override=mappings_override,
                           functional=functional,
                           turn_knobs=turn_knobs,
                           checkpoint_file=checkpoint_file,
                           cal_logger=cal_logger)

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
            slab_struct = SlabGenerator(initial_structure=strt_structure,
                                        miller_index=self.system['hkl'],
                                        min_slab_size=thickness,
                                        min_vacuum_size=vacuum,
                                        lll_reduce=False, center_slab=True,
                                        primitive=False).get_slab()
        slab_struct.sort()
        sd = self.set_sd_flags(slab_struct)
        comment = 'VAC' + str(vacuum) + 'THICK' + str(thickness)
        return Poscar(slab_struct, comment=comment,
                      selective_dynamics=sd)

    @staticmethod
    def set_sd_flags(interface=None, n_layers=2, top=True, bottom=True):
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
        z_coords = slab.frac_coords[:, 2]
        z_lower_bound = None
        z_upper_bound = None
        if bottom:
            z_lower_bound = np.unique(z_coords)[n_layers - 1]
            sd_flags[[i for i, coords in enumerate(slab.frac_coords)
                      if coords[2] <= z_lower_bound]] = np.ones((1, 3))
        if top:
            z_upper_bound = np.unique(z_coords)[-n_layers]
            sd_flags[[i for i, coords in enumerate(slab.frac_coords)
                      if coords[2] >= z_upper_bound]] = np.ones((1, 3))
        return sd_flags.tolist()

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
                 is_matrix=False, Grid_type='A',
                 parent_job_dir='.', job_dir='./Interface',
                 qadapter=None, job_cmd='qsub', wait=True,
                 mappings_override=None, functional="PBE",
                 turn_knobs={'VACUUM': [], 'THICKNESS': []},
                 from_ase=False, checkpoint_file=None,
                 cal_logger=None):
        CalibrateSlab.__init__(self, incar, poscar, potcar, kpoints,
                               system=system, is_matrix=is_matrix,
                               Grid_type=Grid_type,
                               parent_job_dir=parent_job_dir,
                               job_dir=job_dir, qadapter=qadapter,
                               job_cmd=job_cmd, wait=wait,
                               mappings_override=mappings_override,
                               functional=functional,
                               turn_knobs=turn_knobs,
                               from_ase=from_ase,
                               checkpoint_file=checkpoint_file,
                               cal_logger=cal_logger)
        self.interface_setup(turn_knobs=turn_knobs)

    def interface_setup(self, turn_knobs=None):
        if not self.system.get('ligand'):
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
                          ligand=Ligand.from_dict(self.system['ligand']),
                          from_ase=self.from_ase)
        iface.sort()
        sd = self.set_sd_flags(iface, n_layers=2)
        # if there are other paramters that are being varied
        # change the comment accordingly
        comment = self.system['hkl'] + self.system['ligand']['name']
        return Poscar(iface, comment=comment,
                      selective_dynamics=sd)


if __name__ == '__main__':
    # STRUCTURE
    a0 = 3.965
    lvec = [[0.5, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5]]
    lvec = np.array(lvec) * a0
    lattice = Lattice(lvec)
    structure = Structure(lattice, ['Pt'], [[0.0, 0.0, 0.0]],
                          coords_are_cartesian=False)

    # INITIAL VASP INPUT SET
    incarparams = {'System': 'test',
                   'ENCUT': 400,
                   'ISMEAR': 1,
                   'SIGMA': 0.1,
                   'EDIFF': 1E-6}
    incar = Incar(params=incarparams)
    poscar = Poscar(structure, comment='test')
    potcar = Potcar(symbols=poscar.site_symbols, functional='PBE',
                    sym_potcar_map=None)
    kpoints = Kpoints.monkhorst_automatic(kpts=(16, 16, 16),
                                          shift=(0, 0, 0))

    # CALIBRATION INPUT
    system = {'hkl': [1, 1, 1], 'ligand': None}
    turn_knobs = OrderedDict([
        ('SIGMA', [0.025, 0.50]),
        ('POTCAR', [{'Pt': 'Pt'}, {'Pt': 'Pt_pv'}, {'Pt': 'Pt_GW'}]),
        ('IBRION', [1, 2, 3]),
        ('KPOINTS', [k for k in range(20, 40, 10)]),
        ('ENCUT', range(400, 700, 100)),
        ('VACUUM', [10, 12, 15]),
        ('THICKNESS', [11])
    ])
    is_matrix = True
    job_dir = 'Slab'
    job_cmd = ['ls', '-lt']

    # SETUP AND RUN
    cal = CalibrateSlab(incar, poscar, potcar, kpoints,
                        system=system,
                        is_matrix=is_matrix,
                        job_dir=job_dir,
                        turn_knobs=turn_knobs)
    cal.setup()
    cal.run(job_cmd)
