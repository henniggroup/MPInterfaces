from __future__ import division, unicode_literals, print_function

"""
combines instrument, calibrate and interfaces to perform the calibration
and run the actual jobs
"""

import sys
import shutil
import os
import json
import logging

import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints

from mpinterfaces.calibrate import Calibrate, CalibrateMolecule,\
      CalibrateSlab, CalibrateBulk
from mpinterfaces.interface import Interface

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


class Measurement(object):
    """
    Takes in calibrate objects and use that to perform various 
    measuremnt(solvation, binding energy etc) calcuations
    
    sets up and runs static calculation

    Override this class to for custom measuremnts
    """
    def __init__(self, cal_objs, setup_dir='.', parent_job_dir='.',
                 job_dir='./Measurement'):
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
            obj.old_jobs = obj.jobs
            obj.jobs = []
            obj.old_job_dir_list = cal.job_dir_list
            obj.job_dir_list = []
            if isinstance(obj, CalibrateMolecule):
                self.calmol.append(obj)
            elif isinstance(obj, CalibrateSlab):
                self.calslab.append(obj)
            elif isinstance(obj, CalibrateBulk):
                self.calbulk.append(obj)                

    def setup(self):
        """
        setup static jobs for the calibrate objects
        copies CONTCAR to POSCAR
        set NSW = 0
        """
        for cal in self.cal_objs:
            if cal.calc_done:
                cal.incar['NSW'] = 0
                for i, dir in enumerate(cal.old_job_dir_list):
                    job_dir = self.job_dir+os.sep+ \
                        cal.old_jobs[i].name.replace(os.sep, '_').replace('.', '_')+ \
                        os.sep+'STATIC'
                    contcar_file = dir+os.sep+'CONTCAR'            
                    cal.poscar = Poscar.from_file(contcar_file)
                    cal.add_job(job_dir=job_dir)
            else:
                logger.warn('previous calc not done yet or is still running')
                logger.warn('Not setting up the measurement job\n')

    def run(self, job_cmd=None):
        """ run jobs """
        for cal in self.cal_objs:
            if cal.calc_done:
                cal.run()
            elif not cal.isrunning:
                cal.setup()
                cal.run()
            else:
                logger.warn('calibration calc still running')
                logger.warn('try again later')

    def get_energy(self, cal):
        """
        measures the energy of a single cal object
        """
        cal.energies = []
        for job_dir in cal.job_dir_list:
            drone = MPINTVaspDrone(inc_structure=True, 
                                   inc_incar_n_kpoints=False)
            bg =  BorgQueen(drone)
            #bg.parallel_assimilate(rootpath)        
            bg.serial_assimilate(job_dir)
            allentries =  bg.get_data()
            for e in allentries:
                if e:
                    cal.energies.append(e.energy)
                    logger.info(e.energy)

    def make_measurements(self):
        """
        """
        for cal in self.cal_objs:
            self.measure_energy(cal)


class MeasurementSolvation(Measurement):
    """
    Solvation
    """
    def __init__(self, cal_objs, setup_dir='.', parent_job_dir='.', job_dir='./MeasurementSolvation',
                 sol_params={'EB_K':80, 'TAU':0}):
        Measurement.__init__(self, cal_objs=cal_objs, setup_dir=setup_dir, 
                            parent_job_dir=parent_job_dir, job_dir=job_dir):
        self.sol_params = sol_params


    def setup(self):
        """
        setup solvation jobs for the calibrate objects
        copies WAVECAR
        and
        sets the solvation params in the incar file
        works only for cal objects that does only single calc
        """
        for cal in self.cal_objs:
            if cal.calc_done:       
                job_dir = self.job_dir+os.sep+ \
                    cal.old_jobs[0].name.replace(os.sep, '_').replace('.', '_')+ \
                    os.sep + 'SOL'       
                cal.incar['LSOL'] = '.TRUE.'
                for k, v in self.sol_params:
                    cal.incar[k] = v
                if not os.path.exists(job_dir):            
                    os.makedirs(job_dir)
                wavecar_file = cal.old_job_dir_list[0]+os.sep+'WAVECAR'
                shutil.copy(wavecar_file, job_dir+os.sep+'WAVECAR')
                cal.add_job(job_dir=job_dir)
            else:
                logger.warn('previous calc in the dir, '+cal.job_dir+'not done yet or is still running')
                logger.warn('Not setting up the measurement job\n' )

    def make_measurements(self):
        """
        call this after setup and run
        get solvation energies
        """
        pass


class MeasurementInterface(Measurement):
    """
    Interface
    """
    def __init__(self, cal_objs, setup_dir='.', parent_job_dir='.', job_dir='./MeasurementSolvation'):
        Measurement.__init__(self, cal_objs=cal_objs, setup_dir=setup_dir, 
                            parent_job_dir=parent_job_dir, job_dir=job_dir):

    def setup(self):
        """
        setup static jobs for the calibrate objects
        copies CONTCAR to POSCAR
        set NSW = 0
        write system.json file for database
        """
        d = {}
        for cal in self.cal_objs:
            if cal.calc_done:
                cal.incar['NSW'] = 0
                for i, dir in enumerate(cal.old_job_dir_list):
                    job_dir = self.job_dir+os.sep+ \
                        cal.old_jobs[i].name.replace(os.sep, '_').replace('.', '_')+ \
                        os.sep+'STATIC'
                    contcar_file = dir+os.sep+'CONTCAR'            
                    cal.poscar = Poscar.from_file(contcar_file)
                    if isinstance(cal, CalibrateSlab) or isinstance(obj, CalibrateInterface):
                        d['hkl'] = cal.hkl
                    if isinstance(obj, CalibrateInterface):                    
                        d['ligand'] = cal.ligand.composition.formula
                    os.makedirs(job_dir)
                    if d:
                        with open(job_dir+os.sep+'system.json') as f:
                            json.dump(d, f)
                    cal.add_job(job_dir=job_dir)
            else:
                logger.warn('previous calc not done yet or is still running')
                logger.warn('Not setting up the measurement job\n')


    def make_measurements(self):
        """
        call this after setup and run
        """
        Binding_energy= self.calculate_binding_energy()
        logger.info(Binding_energy)

    def calculate_binding_energy(self): 
        """
        calculates the binding energies as Binding Energy = Interface - (Slab + Ligand) 
        in cal_objs list , 0th is treated as interface, 1st is treated as slab, 2nd 
        is treated as ligand
        """
        self.get_energy(self.cal_objs)
        Interface_energy= self.cal_objs[0].energies[0]
        Slab_energy= self.cal_objs[1].energies[1]
        Ligand_energy= self.cal_objs[2].energies[2]
        #testing purpose 
        if Ligand_energy == None or Slab_energy == None or Interface_energy == None:
                return "Binding Energy not ready"
        else:
                Binding_energy= Interface_energy - (Slab_energy + n_ligands*Ligand_energy)
                return Binding_energy


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
    measure.setup()
    measure.run()
    #set the measurement jobs
    #for cal in cal_objs:                    
    #    measure.setup(cal)
    #    #measure.setup_solvation_job(cal)
    #measure.run()

