# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

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

from six.moves import range

import sys
import os
import re
import datetime
from itertools import product
from collections import OrderedDict

import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.inputs import Potcar, Kpoints
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.symmetry.bandstructure import HighSymmKpath

from custodian.custodian import Custodian

from monty.json import MSONable, MontyEncoder
from monty.serialization import dumpfn

from mpinterfaces.instrument import MPINTVaspInputSet, MPINTVaspJob
from mpinterfaces.interface import Interface, Ligand
from mpinterfaces.utils import get_ase_slab, get_magmom_string, get_magmom_afm, \
    get_magmom_mae, print_exception
from mpinterfaces.mat2d.electronic_structure import get_2D_hse_kpoints,\
    get_2D_incar_hse_prep, get_2D_incar_hse
from mpinterfaces.default_logger import get_default_logger

__author__ = "Kiran Mathew, Joshua J. Gabriel"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


logger = get_default_logger(__name__)


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
                 database=None, magnetism=None, mag_init=None, reuse=None,
                 reuse_override=None, reuse_incar=None, solvation=None,
                 turn_knobs=OrderedDict([('ENCUT', []),
                                         ('KPOINTS', [])]),
                 checkpoint_file=None, finer_kpoint=None, cal_logger=None):
        """
        Calibrate constructor

        Args:
            incar (Incar object): input INCAR
            poscar (Poscar object): input POSCAR
            potcar (Potcar object): input POTCAR
            kpoints (Kpoints object): input KPOINTS
            system: system info as a dictionary,
                slab or interface example:
                system={'hkl':[1,1,1], 'ligand':None},
            is_matrix (bool): whether the jobs are dependent on each
                other
            Grid_type (str): kpoints grid_type
            parent_job_dir (str): the directory from which all the
                jobs are launched
            job_dir (str): job directory created for each job in the
                parent_job_dir
            qadapter (?): adapter for the batch system
            job_cmd (str): command to be used for submitting the job. If
                qadapter is specified then job_cmd is ignored
            wait (bool): whther to wait for the job to finish. If the job is
                being submitted to the queue then there is no need for
                waiting
            turn_knobs (dict): an ordered dictionary of parmaters and the
                corresponding values
            mappings_override (dict): override symbol mapping in potcar
                               eg:- {'S':'S_sv'}
            functional (str): exchange-correlation functional
            database (str): A work in progress, will be a database_name.yaml
                            file for defaults specific to a database workflow
                            that will have defaults for
                            INCAR: cutoff, convergence for relaxation and
                                   continuation jobs
                            KPOINTS: for relaxation, band structure jobs
                            POTCAR: database specific
                            For now defaults to None, if set to 'twod'
                            activates twod set of directives
            reuse (list or bool): list of filenames for reuse
                          Eg: ['CHGCAR', 'WAVECAR']
                          'CONTCAR' is copied by default and if found empty
                          warning is issued. Use the following flag for override
                          only if you know what you are doing
                          'True' for just copying the CONTCAR file
            reuse_override (bool): whether to override the missing CONTCAR for a
                          reuse calc
            magnetism (str): specifies magnetism calculation to be used
                           implemented are 'AntiFerroMagnetism' and
                           'Magntic Anisotropy Energy'
            solvation (bool): whether to activate a solvation job, sets LSOL=True
                           for now

        Calibrate jobs represent the engine configuration of mpinterfaces,
        where the fuel (input file sources) and driving method (kind of calculation)
        are decided . The Engine itself is instrument.py which creates the input set
        configured in Calibrate.

        Current fueling methods:
          1. simplest test case involving a single job:
               - specify the incar, kpoints, poscar, potcar (aka the VASP 4)
                 explicitly as pymatgen objects
               - turn_knobs = {} , is_matrix = False
          2. test case for calibration of parameters:
               - specify an initial configuration for the VASP 4
               - specify parameters to calibrate via turn_knobs,
                 set is_matrix = True only if number of parameters > 1
          3. Database production case: (possibly most used)
               - specify initial configuration for the VASP 4 based on
                 a database.yaml
               - specify an input.yaml that details the workflow

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
        self.database = database
        self.magnetism = magnetism
        self.mag_init = mag_init
        self.solvation = solvation
        self.reuse = reuse
        self.reuse_incar = reuse_incar
        self.reuse_override = reuse_override
        self.reuse_paths = None  # list object communicated to instrument
        self.finer_kpoint = finer_kpoint
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
                    self.setup_potcar_jobs(mappings=v, functional_list=None)
                elif k == 'POTCAR_functional' and v:
                    self.setup_potcar_jobs(mappings=None, functional_list=v)
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
        n_items = len(list(self.turn_knobs.items()))
        keys = list(self.turn_knobs.keys())
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
        #### Testing ####
        # Orig
        # job_dir = self.job_dir + os.sep + self.key_to_name(keys[i])
        job_dir = '__'.join(
            [self.job_dir.split('/')[-1], self.key_to_name(keys[i])])
        #### Testing ####
        if i == n - 1 and i != 0:
            for val in self.turn_knobs[keys[i]]:
                ##
                ## self.job_dir = job_dir + os.sep + self.val_to_name(val)
                self.job_dir = '__'.join([job_dir, self.val_to_name(val)])
                self.logger.info(
                    'setting jobs in the directory: ' + self.job_dir)
                self._setup(turn_knobs=dict([(keys[i], [val])]))
                self.add_job(name=job_dir, job_dir=self.job_dir)
        else:
            for val in self.turn_knobs[keys[i]]:
                ##
                ## self.job_dir = job_dir + os.sep + self.val_to_name(val)
                self.job_dir = '__'.join([job_dir, self.val_to_name(val)])
                self.logger.info(
                    'setting jobs in the directory: ' + self.job_dir)
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
        if isinstance(val, float):
            return re.sub('\.', '_', str(val))
        elif isinstance(val, list):
            return self.kpoint_to_name(val, 'M')
        elif isinstance(val, dict):
            return self.potcar_to_name(mapping=val)
        elif isinstance(val, Poscar):
            name = str(val.structure.composition.reduced_formula) \
                + '_' + str(int(val.structure.lattice.volume)) \
                + '_' + ''.join((val.comment).split())
            return name.replace('\\', '_').replace('(', '_').replace(')', '_')
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
            kpoint = [str(k).replace('.', '_') for k in kpoint]
            return 'x'.join(kpoint)
        else:
            return str(kpoint)

    def potcar_to_name(self, mapping=None, functional=None):
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
            return '_'.join(l)
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
                if sym in list(self.mappings_override.keys()):
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

    def set_kpoints(self, kpoint=None, poscar=None, ibzkpth=None):
        """
        set the kpoint
        """
        # useful to check if a poscar is supplied from setup_poscar_jobs (most often the case)
        # or this is a single poscar use case
        if not poscar:
            poscar = self.poscar

        # splitting into two if elif branches means fewer if statements to check on
        # a run

        # Most general method of setting the k-points for
        # different grid types
        # NOTE: requires that at least one k-points value be passed
        # as a turn - knobs list value
        # this is not true for values that may be caculated out of
        # a database
        # use this part only if this is a non-database run for example
        # for k-points calibration

        if not self.database:

            if self.Grid_type == 'M':
                self.kpoints = Kpoints.monkhorst_automatic(kpts=kpoint)
            elif self.Grid_type == 'A':
                self.kpoints = Kpoints.automatic(subdivisions=kpoint)
            elif self.Grid_type == 'G':
                self.kpoints = Kpoints.gamma_automatic(kpts=kpoint)
            elif self.Grid_type == '3D_vol':
                self.kpoints = Kpoints.automatic_density_by_vol(structure=poscar.structure,
                                                                kppvol=kpoint)
            elif self.Grid_type == 'bulk_bands_pbe':
                self.kpoints = Kpoints.automatic_linemode(divisions=kpoint,
                                                          ibz=HighSymmKpath(
                                                              poscar.structure))

            elif self.Grid_type == 'D':
                self.kpoints = Kpoints.automatic_density(structure=poscar.structure,kppa=kpoint)

            elif self.Grid_type == 'Finer_G_Mesh':
                # kpoint is the scaling factor and self.kpoints is the old kpoint mesh
                self.logger.info('Setting Finer G Mesh for {0} by scale {1}'.format(kpoint, self.finer_kpoint))
                self.kpoints = Kpoints.gamma_automatic(kpts = \
                   [i * self.finer_kpoint for i in kpoint])
                self.logger.info('Finished scaling operation of k-mesh')

        # applicable for database runs
        # future constructs or settinsg can be activated via a yaml file
        # database yaml file or better still the input deck from its speification
        # decides what combination of input calibrate constructor settings to use
        # one of them being the grid_type tag

        elif self.database == 'twod':

            # set of kpoints settings according to the 2D database profile
            # the actual settings of k-points density
            # will in future come from any database input file set

            if self.Grid_type == 'hse_bands_2D_prep':
                kpoint_dict = Kpoints.automatic_gamma_density(poscar.structure,
                                                              200).as_dict()
                kpoint_dict['kpoints'][0][2] = 1  # remove z kpoints
                self.kpoints = Kpoints.from_dict(kpoint_dict)

            elif self.Grid_type == 'hse_bands_2D':
                # can at most return the path to the correct kpoints file
                # needs kpoints to be written out in instrument in a different way
                # not using the Kpoints object
                self.kpoints = get_2D_hse_kpoints(poscar.structure, ibzkpth)

            elif self.Grid_type == 'bands_2D':
                kpoint_dict = Kpoints.automatic_linemode(divisions=20,
                                                         ibz=HighSymmKpath(poscar.structure)).as_dict()
                self.kpoints = Kpoints.from_dict(kpoint_dict)

            elif self.Grid_type == 'relax_2D':
                # general relaxation settings for 2D
                kpoint_dict = Kpoints.automatic_gamma_density(poscar.structure,
                                                              1000).as_dict()
                kpoint_dict['kpoints'][0][2] = 1
                self.kpoints = Kpoints.from_dict(kpoint_dict)

            elif self.Grid_type == 'relax_3D':
                # general relaxation settings for 3D
                kpoint_dict = Kpoints.automatic_gamma_density(
                    poscar.structure, 1000)
                self.kpoints = Kpoints.from_dict(kpoint_dict)

    def setup_incar_jobs(self, param, val_list):
        """
        set up incar jobs,, calls set_incar to set the value to param

        Args:
            param: Name of INCAR parameter
            val_list: List of values to vary for the param
        """

        if val_list != ['2D_default']:
            for val in val_list:
                self.logger.info('setting INCAR parameter ' + param + ' = '
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
                    job_dir = self.job_dir + os.sep + self.key_to_name(
                        'KPOINTS') \
                        + os.sep + self.kpoint_to_name(kpoint,
                                                       self.Grid_type)
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
            for pos in poscar_list:
                # if it is a twod_database run or any general standard database run,
                # the incar, kpoints and potcar follow a standard input set
                # which will be activated by the twod_database tag set to true
                # NOTE: this implementation means that the first turn_knobs tag
                # needs to be the poscar objects list
                # the database tag will be set to the name of the yaml file with the
                # standard input deck definition for that database
                # this incar dict provided as the init can be general format
                # based on the chosen functional, cutoff
                # so self.incar is a vdW incar for re-relaxation in vdW, gga for every
                # other calculation or LDA+U for LSDA+U calculations
                incar_dict = self.incar.as_dict()
                if self.reuse:
                    # if this is a true list minimally, ['CONTCAR']
                    # it is to be ensured that the poscar list is a
                    # list of paths as opposed to list of poscar objects by the turn knobs
                    # values
                    # Here pos is the path and r in each self.reuse is the name of the file(s)
                    # to be reused
                    # in a reuse calculation the following are possible:
                    # update incar (Solvation calculations) or reset incar (HSE calculations)
                    # reset kpoints file with IBZKPT
                    # copy a CHGCAR or WAVECAR or both perhaps
                    try:
                        # first setup of POSCAR initial, INCAR, KPOINTS
                        poscar = Poscar.from_file(pos + os.sep + 'CONTCAR')
                        self.logger.info('Read previous relaxed CONTCAR file from {}'.
                                         format(pos))
                        # check if it is KPOINTS altering job like HSE
                        if self.Grid_type == 'hse_bands_2D_prep':
                            # HSE prep calcualtions
                            # reset the INCAR file with a magmom only if exists
                            try:
                                incar_dict = {
                                    'MAGMOM': get_magmom_string(poscar)}
                            except:
                                incar_dict = {}
                            incar_dict = get_2D_incar_hse_prep(incar_dict)
                            self.set_kpoints(poscar=poscar)
                            self.logger.info(
                                'updated input set for HSE 2D prep calcaultion')

                        elif self.Grid_type == 'hse_bands_2D':
                            # HSE calculation
                            # reset the incar and kpoints file builds
                            # on the preceding calculations (prep calculation)
                            # IBZKPT
                            try:
                                incar_dict = {
                                    'MAGMOM': get_magmom_string(poscar)}
                            except:
                                incar_dict = {}
                            incar_dict = get_2D_incar_hse(incar_dict)
                            self.set_kpoints(poscar=poscar,
                                             ibzkpth=pos + os.sep + 'IBZKPT')
                            self.logger.info('updated input set for HSE calcaultion\
                                             using IBZKPT from {0}'.format(pos + os.sep + 'IBZKPT'))

                        elif self.Grid_type == 'hse_bands':
                            # general HSE bands
                            pass

                        elif self.Grid_type == 'Finer_G_Mesh':
                            self.logger.info('updating to Finer G Mesh')
                            kpoint = Kpoints.from_file(pos+os.sep+'KPOINTS')
                            self.set_kpoints(kpoint=kpoint.kpts[0])

                        else:
                            # use the same kpoints file and build from the old
                            # incar
                            self.kpoints = Kpoints.from_file(
                                pos + os.sep + 'KPOINTS')
                            # decide on how to use incar, use same one or
                            # update or afresh
                            if self.reuse_incar == 'old':
                                incar_dict = Incar.from_file(
                                    pos + os.sep + 'INCAR').as_dict()
                            elif self.reuse_incar == 'update':
                                # way to go for cutoff updates, convergence, etc.
                                # but retain the old functional
                                incar_dict.update(Incar.from_file(pos + os.sep + 'INCAR').
                                                  as_dict())
                            else:
                                # use a fresh incar as specified by the init
                                # way to go for example for LDAU or other
                                # major removals done to INCAR
                                # but always retain the MAGMOM if present
                                old_incar_dict = Incar.from_file(
                                    pos + os.sep + 'INCAR').as_dict()
                                if 'MAGMOM' in old_incar_dict.keys():
                                    incar_dict['MAGMOM'] = old_incar_dict[
                                        'MAGMOM']
                                else:
                                    incar_dict = incar_dict

                        if isinstance(self.reuse, list):
                            reuse_paths = [
                                pos + os.sep + r for r in self.reuse]
                            self.reuse_paths = reuse_paths
                        # Magnetism use cases, updates to be made to the INCAR (MAE)
                        # and poscar (AFM)
                        # MAE and AFM

                        if self.magnetism == 'MAE':

                            # remove vdW tags for MAE calculations
                            vdW_tags = ('GGA', 'AGGAC', 'LUSE_VDW',
                                        'PARAM1', 'PARAM2')
                            for key in vdW_tags:
                                if key in incar_dict:
                                    del incar_dict[key]

                            self.logger.info(
                                'updating input set for MAE calculation')
                            self.mag_init = Outcar(
                                pos + os.sep + 'OUTCAR').total_mag
                            nbands = 2 * \
                                Vasprun(pos + os.sep +
                                        'vasprun.xml').parameters['NBANDS']
                            # u_value = Vasprun(pos+os.sep+'vasprun.xml').incar['LDAUU']
                            # u_value = 4.0
                            self.logger.info(
                                "updating mag mom with value {0}".format(self.mag_init))
                            self.logger.info(
                                "updating NBANDS with {0}".format(nbands))
                            incar_dict.update({'NBANDS': nbands,
                                               'LSORBIT': True,
                                               'EDIFF': 1e-08,
                                               'ICHARG': 11,
                                               'LMAXMIX': 4,
                                               'LCHARG': False,
                                               'ISYM': 0,
                                               'NSW': 0,
                                               'ISPIN': 2,
                                               'IBRION': -1,
                                               'LORBIT': 11,
                                               'MAGMOM': get_magmom_mae(poscar, self.mag_init)
                                               })
                            # incar_dict.update({'LDAUU': u_value})

                        elif self.magnetism == 'AFM':
                            self.logger.info(
                                'updating INCAR and POSCAR for AFM calculation')
                            afm, poscar = get_magmom_afm(poscar, self.database)
                            incar_dict.update({'MAGMOM': afm})
                    except:
                        # check what to do if the previous calculation being reused is not
                        # actuall done .. system exit or adopt a user override
                        # with POSCAR
                        self.logger.warn(
                            'Empty relaxed CONTCAR file .. Probably job not done')
                        if not self.reuse_override:
                            self.logger.warn(
                                'You can set reuse_override to continue with POSCAR file, exiting now ..')
                            sys.exit(0)
                        else:
                            self.logger.info('Using old Poscar for rerun')
                            poscar = Poscar.from_file(pos + os.sep + 'POSCAR')

                # case for non - reuse
                else:
                    poscar = pos
                    # temporary: magnetism only set if twod flag is activated
                    if self.database == 'twod':
                        incar_dict.update(
                            {'MAGMOM': get_magmom_string(poscar)})
                        self.set_kpoints(poscar=poscar)
                    self.incar = Incar.from_dict(incar_dict)

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
                                self.qadapter, vis_logger=self.logger,
                                reuse_path=self.reuse_paths)
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
        z_coords = interface.frac_coords[:, 2]
        z_coords_slab = slab.frac_coords[:, 2]
        z_lower_bound = None
        z_upper_bound = None
        if bottom:
            z_lower_bound = np.unique(z_coords_slab)[n_layers - 1]
            sd_flags[np.where(z_coords <= z_lower_bound)] = np.ones((1, 3))
        if top:
            z_upper_bound = np.unique(z_coords_slab)[-n_layers]
            sd_flags[np.where(z_coords >= z_upper_bound)] = np.ones((1, 3))
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
        iface = sorted(Interface(structure,
                                 hkl=self.system['hkl'],
                                 ligand=Ligand.from_dict(
                                     self.system['ligand']),
                                 from_ase=self.from_ase))
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
        ('ENCUT', list(range(400, 700, 100))),
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
