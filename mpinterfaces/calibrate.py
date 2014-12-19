"""
Calibration module

"""

import sys
import os, shutil
import shlex, subprocess
import operator
from collections import Counter
import re
import time
import datetime
from pprint import pprint
import logging
import numpy as np
from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
#from custodian.vasp.handlers import VaspErrorHandler, FrozenJobErrorHandler, MeshSymmetryErrorHandler, NonConvergingErrorHandler
from custodian.custodian import Custodian, gzip_dir
from custodian.vasp.interpreter import VaspModder
from pymatgen.apps.borg.queen import BorgQueen
from mpinterfaces.instrument import MPINTVaspInputSet, MPINTVaspJob
from mpinterfaces.data_processor import MPINTVaspDrone


class Calibrate(object):
    
    """
        
    The base class for creating vasp work flows for
    calibrating the input parameters for different systems
    
    
    """

    def __init__(self, incar, poscar, potcar, kpoints,
                 setup_dir='.', parent_job_dir='.',job_dir='./Job',
                 qadapter=None, job_cmd='qsub', turn_knobs={'ENCUT':[],'KPOINTS':[]} ):
        """
        setup_dir = directory from where the setup files are copied from
        parent_job_dir = the directory from which the script is run
        job_dir = name of the job directory relative to the parent directory
        """
        self.parent_job_dir = os.path.abspath(parent_job_dir)
        self.setup_dir = os.path.abspath(setup_dir)
        self.job_dir = job_dir
        self.incar = incar
        self.poscar = poscar
        self.potcar =potcar
        self.kpoints = kpoints
        self.qadapter = qadapter
        self.vis = []
        self.jobs = []
        #example:- handlers = [VaspErrorHandler(), FrozenJobErrorHandler(),
        #MeshSymmetryErrorHandler(), NonConvergingErrorHandler()]        
        self.handlers = [ ]
        self.job_cmd = job_cmd
        self.n_atoms = 0
        self.turn_knobs = turn_knobs
        self.response_to_knobs = {}
        for k, v in turn_knobs.items():
            self.response_to_knobs[k] = {}
            


    def setup(self):
        for k, v in self.turn_knobs.items():
            if k == 'KPOINTS' and v:
                self.setup_kpoints_jobs(kpoints_list = v)
            elif v:
                self.setup_incar_jobs(k, v)                
                

    def setup_incar_jobs(self, param, val_list):
        for val in val_list:
            print 'setting '+param+' = ', val
            job_dir  = self.job_dir+ os.sep + param + os.sep + str(val)
            self.incar[param] = val
            self.add_job(name=param+str(val), job_dir=job_dir)
                
            
    def setup_kpoints_jobs(self, Grid_type = 'M', kpoints_list = None, conv_step = 1):
        """
        set the jobs for kpoint convergence
        
        """
        if Grid_type == 'M':
            #local list convergence_list , convert from tuple
            #because constructor takes tuple as argument
            if kpoints_list:
                conv_list = kpoints_list 
                start = conv_list[0]
                end = conv_list[1]
                if (conv_step):
                    for x in range(1+start[0], end[0], conv_step):
                        conv_list.append([x, x, x])
                    for kpoint in conv_list:
                        self.kpoints = Kpoints.monkhorst_automatic(kpts = kpoint)
                        name = str(kpoint[0]) + 'x' + str(kpoint[1]) + 'x' + str(kpoint[2])
                        print 'KPOINTS = ', name
                        job_dir = self.job_dir +os.sep+ 'KPOINTS' + os.sep + name
                        self.add_job(name=name, job_dir=job_dir) 
            else:
                print 'kpoints_list not provided'


        
    def add_job(self, name='noname', job_dir='.'):
        vis = MPINTVaspInputSet(name, self.incar, self.poscar,
                                self.potcar, self.kpoints, self.qadapter)
        #the job command can be overrridden in the run method
        job = MPINTVaspJob(self.job_cmd.split(), final = True, setup_dir=self.setup_dir,
                               parent_job_dir=self.parent_job_dir, job_dir=job_dir,
                               vis=vis, auto_npar=False, auto_gamma=False)
        self.jobs.append(job)
                

    
    def run(self, job_cmd=None):
        """
        run the vasp jobs through custodian

        """
        #if the job list is empty, run a single job with the provided input set
        if not self.jobs :
            self.add_job(name='single job', job_dir=self.job_dir)             
        
        #override the job_cmd if provided
        for j in self.jobs:
            if job_cmd is not None:            
                j.job_cmd = job_cmd
            else:
                j.job_cmd = self.job_cmd

        c_params = {'jobs': [j.as_dict() for j in self.jobs],
                'handlers': [h.as_dict() for h in self.handlers], 'max_errors': 5}
        c = Custodian(self.handlers, self.jobs, max_errors=5)
        c.run()

    def get_knob_responses(self):
        drone = MPINTVaspDrone(inc_structure=True, inc_incar_n_kpoints=True) 
        bg =  BorgQueen(drone)
        for k, v in self.response_to_knobs.items():
            rootpath = self.job_dir+ os.sep + k
            print 'rootpath', rootpath
            #bg.parallel_assimilate(rootpath)        
            bg.serial_assimilate(rootpath)
            allentries =  bg.get_data()
            print allentries
            for e in allentries:
                if e:
                    self.n_atoms = len(e.structure)
                    print 'n_atoms', self.n_atoms
                    if k == 'KPOINTS':
                        self.response_to_knobs[k][str(e.kpoints.kpts)] = e.energy
                    else:
                        self.response_to_knobs[k][str(e.incar[k])] = e.energy
        


    def enforce_cutoff(self, input_list, delta_e_peratom=0.001):
        """
        energy difference of 1meV per atom
        """
        matching_list = []
        for i, e in enumerate(input_list):
            if i < len(input_list)-1:
                print i, input_list[i+1], e
                print np.abs(input_list[i+1][1] - e[1])
                if np.abs(input_list[i+1][1] - e[1])/self.n_atoms <= delta_e_peratom:
                    matching_list.append(input_list[i+1][0])
        if matching_list:
            print matching_list
            matching_kpt_list = []
            if '[[' in matching_list[0]:
                for ml in matching_list:
                    if '[[' in ml:
                        m = re.search(r"\[\[(\d+)\,(\d+)\,(\d+)\]\]", ml)
                        matching_kpt_list.append( [ int(m.group(1)), int(m.group(2)), int(m.group(3))])
                return matching_kpt_list
            else:
                return [float(encut) for encut in matching_list]
        else:
            return []
        

            
    def optimum_params(self, allentries, en_mc, kp_mc):
        """
        input: all enetires, values of encut and kpoints used for kpoints and encut studies rexpectively
        sets the dictionaries of kpoints and energies and  encut and energies
        returns: optimum kpoints and encut
        """
        #dict of kpoints and energies
        kpt={}
        #dict of encut and energies
        encut = {}
        for e in allentries:
            if e:
                if str(e.incar['ENCUT']) == en_mc:
                    kpt[str(e.kpoints.kpts)] = e.energy
                if str(str(e.kpoints.kpts)) == kp_mc:
                    encut[str(e.incar['ENCUT'])] = e.energy
        #order the keys(encut or kpoint)from large value of the energy  to small value
        sorted_encut = sorted(encut.items(), key=operator.itemgetter(1), reverse=True)
        sorted_kpt = sorted(kpt.items(), key=operator.itemgetter(1), reverse=True)
        print encut,'\n', sorted_encut
        print kpt, '\n', sorted_kpt
        #get the list of encut and kpoints that satisfy the delate criterion
        #mind: default deltae = 0.01eV
        matching_encut = self.enforce_cutoff(sorted_encut)
        matching_kpt = self.enforce_cutoff(sorted_kpt)
        opt_encut = None        
        opt_kpt = None
        #of the possible encuts and kpoints, pick the optimum one
        #i.e for encut, the lowest value and for kpoints the one that corresponds
        #to the lowest number of kpoints
        if matching_encut:
            opt_encut = np.min(np.array(matching_encut))
        else:
            print 'no ENCUT met the convergence criterion'
        if matching_kpt:
            nkpt = matching_kpt[0][0] * matching_kpt[0][1] * matching_kpt[0][2]
            for i, val in enumerate(matching_kpt):
                if i < len(matching_kpt)-1:
                    nkpt1 = val[i+1][0] * val[i+1][1] * val[i+1][2]
                    if nkpt1<nktp:
                        opt_kpt = val
        else:
            print 'no KPOINTS met the convergence criterion'
        #opt_encut: list of floats
        #opt_kpt: list of list of integers
        return opt_encut, opt_kpt
                                

        
    def knob_settings(self, rootpath=None):
        """
        go through the parent dir and get all encut, kpoints and energies
        also vac spacing and slab thinckness for slab calulations
        these values willl be used to do the actual interface measurements
        use Vasprun class to get the afore mentioned values from the xml files
        should not proceed if the calculations are not done
        should update the incar, poscar, potcar, kpoints objects according to
        the knob_settings
        
        
        """
        if rootpath is None:
            rootpath = self.job_dir
        print 'xx', self.job_dir
        drone = MPINTVaspDrone(inc_structure=True, inc_incar_n_kpoints=True) #VaspToComputedEntryDrone()#
        bg =  BorgQueen(drone)
        #bg.parallel_assimilate(rootpath)        
        bg.serial_assimilate(rootpath)
        allentries =  bg.get_data()
        alldata = []
        for e in allentries:
            if e:
                self.n_atoms = len(e.structure)
                alldata.append(str(e.incar['ENCUT']))
                alldata.append(str(e.kpoints.kpts))
        #get the 2 most common items in alldata
        enkp_mc =  Counter(alldata).most_common(2)
        kp_mc = None
        en_mc = None
        #if the most common item is kpoints then for the encut convergence study, that value of kpoint was used
        #else the other way around
        if '[[' in enkp_mc[0][0]:
            kp_mc = enkp_mc[0][0]
            en_mc = enkp_mc[1][0]
        else:
            kp_mc = enkp_mc[1][0]
            en_mc = enkp_mc[0][0]
        #opt_encut: list of floats
        #opt_kpt: list of list of integers
        opt_encut, opt_kpt = self.optimum_params(allentries, en_mc, kp_mc)
        


        
class CalibrateMolecule(Calibrate):
    
    """
    
    Calibrate paramters for Molecule calculations
    
    """

    def __init__(self, incar, poscar, potcar, kpoints,
                 setup_dir='.', parent_job_dir='.', job_dir='./Molecule',
                qadapter=None, job_cmd='qsub', turn_knobs={'ENCUT':[],'KPOINTS':[]} ):
        
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                           setup_dir=setup_dir, parent_job_dir=parent_job_dir,
                           job_dir=job_dir, qadapter=qadapter, job_cmd=job_cmd,
                           turn_knobs = turn_knobs)


        
    def setup_kpoints_jobs(self, Grid_type = 'M', kpoints_list = None, conv_step = 1):
        print "Its a molecule ! no need for kpoint convergence"
        return

        
                

class CalibrateBulk(Calibrate):
    
    """
    
    Calibrate paramters for Bulk calculations
    
    """
    
    def __init__(self, incar, poscar, potcar, kpoints,
                  setup_dir='.', parent_job_dir='.', job_dir='./Bulk',
                qadapter=None, job_cmd='qsub', turn_knobs={'ENCUT':[],'KPOINTS':[]}): 
            
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                            setup_dir=setup_dir, parent_job_dir=parent_job_dir,
                             job_dir=job_dir, qadapter=qadapter, job_cmd=job_cmd,
                             turn_knobs = turn_knobs )

        

                    
class CalibrateSlab(Calibrate):
    
    """
    
    Calibrate paramters for Slab calculations
    
    """
    
    def __init__(self, incar, poscar, potcar, kpoints,
                 setup_dir='.', parent_job_dir='.', job_dir='./Slab',
                qadapter=None, job_cmd='qsub', turn_knobs={'ENCUT':[],'KPOINTS':[]}):
        
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                            setup_dir=setup_dir, parent_job_dir=parent_job_dir,
                            job_dir=job_dir, qadapter=qadapter, job_cmd=job_cmd,
                            turn_knobs = turn_knobs)

    def setup_kpoints_jobs(self, Grid_type = 'M', kpoints_list = None, conv_step = 1):
        if Grid_type == 'M':
            #local list convergence_list , convert from tuple
            #because constructor takes tuple as argument
            if kpoints_list:
                conv_list = kpoints_list 
                start = conv_list[0]
                end = conv_list[1]
                if (conv_step):
                    for x in range(1+start[0], end[0], conv_step):
                        conv_list.append([x, x, 1])                         
                    for kpoint in conv_list:
                        self.kpoints = Kpoints.monkhorst_automatic(kpts = kpoint)
                        name = str(kpoint[0]) + 'x' + str(kpoint[1]) + 'x' + str(kpoint[2])
                        print 'KPOINTS = ', name
                        job_dir = self.job_dir +os.sep+ 'KPOINTS' + os.sep + name
                        self.add_job(name=name, job_dir=job_dir)
            else:
                print 'kpoints_list not provided'
        
	            


#test
if __name__ == '__main__':

    system = 'Pt bulk'
    atoms = ['Pt']
    
    a0 = 3.965
    lvec = [ [0.5, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5] ]
    lvec = np.array(lvec) * a0
    lattice = Lattice(lvec)#.from_parameters(3.866, 3.866, 3.866, 60, 60, 60)
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
    atoms = poscar.site_symbols
    potcar = Potcar(symbols=atoms, functional='PBE', sym_potcar_map=None)
    kpoints = Kpoints.monkhorst_automatic(kpts=(16, 16, 16), shift=(0, 0, 0))#{'grid_density': 1000} #

#    calmol = CalibrateMolecule(incar, poscar, potcar, kpoints)
    calbulk = CalibrateBulk(incar, poscar, potcar, kpoints, job_dir='./Bulk_test',
                            turn_knobs = {'ENCUT':range(400,800,100),
                                          'KPOINTS':[ [7, 7, 7], [11, 11, 11] ] } )    
    calbulk.setup()
    #the job_cmd can passed to the run
    #['qsub','job_script']
    #calbulk.run(['ls','-lt'])
    
    #get all data in all the directories in the provided rootfolder, here 1/
    calbulk.get_knob_responses()
    print calbulk.response_to_knobs
    #test enforce_cutoff
    #inp_list = [ ['[[2,2,4]]', 10], ['[[2,2,5]]', 9.9], ['[[2,2,6]]', 9.895], ['[[2,2,7]]', 9.888], ['[[2,2,8]]', 9.879],]
    #print calbulk.enforce_cutoff(inp_list, delta_e=0.01)




#note: write the default yaml to the directory where the jobs are run. This is useful later on for comparing different job runs in that direcctory
