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
                 setup_dir='.', parent_job_dir='.',job_dir='.'):
        """
        setup_dir = directory from where the setup files are copied from
        eg:- the submit script for the queue system
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
        self.vis = []
        self.jobs = []
        #example:- handlers = [VaspErrorHandler(), FrozenJobErrorHandler(),
        #MeshSymmetryErrorHandler(), NonConvergingErrorHandler()]        
        self.handlers = [ ] 


        
    def encut_cnvg(self, encut_list):
        """
        create ENCUT convergence workflow
        """
        for encut in encut_list:
            print 'ENCUT = ', encut
            job_dir  = self.job_dir+ os.sep + 'ENCUT' + os.sep + str(encut)
            self.incar['ENCUT'] = encut
            vis = MPINTVaspInputSet('encut_'+str(encut), self.incar,
                                    self.poscar, self.potcar, self.kpoints)
            #the job command can be overrridden in the run method
            job = MPINTVaspJob(["pwd"], final = True, setup_dir=self.setup_dir,
                               parent_job_dir=self.parent_job_dir, job_dir=job_dir,
                               vis=vis, auto_npar=False, auto_gamma=False)
            self.jobs.append(job)


            
    def kpoints_cnvg(self, Grid_type = 'M', kpoints_list = None, conv_step = 1):
        """
        set the jobs for kpoint convergence
        
        Calls Kpoints constructors according to Grid_type:
        Monkhorst Pack Automatic is default,
        
        G for Gamma centered automatic; ONLY MP method implemented now others can be added 
        
        kpoints_list describes the start and end of kpoints set
        eg: user can pass [[6, 6, 6], [10, 10, 10]] and conv_step for a
        convergence to be done for 6x6x6 to 10x10x10,
        
        is_slab can be switched to True constraints z to default to 1,
        that is 6x6x1 to 10x10x1
        
        first defines the list of kpoints according to the user input,
        user needs to give only the start and end kpoint,
        whether it is for a slab or simple bulk
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
                        print 'KPOINTmesh = ', name
                        job_dir = self.job_dir +os.sep+ 'KPOINTS' + os.sep + name
                        vis = MPINTVaspInputSet(name, self.incar, self.poscar, self.potcar,
                                                self.kpoints)
                        job = MPINTVaspJob(["pwd"], final = True, setup_dir=self.setup_dir,
                                parent_job_dir=self.parent_job_dir, job_dir=job_dir,
                                vis=vis, auto_npar=False, auto_gamma=False)
                        
                        self.jobs.append(job)
            else:
                print 'kpoints_list not provided'


    
    def run(self, job_cmd=None):
        """
        run the vasp jobs through custodian
        set up the custodian that we want to run : just a single run
        just one vasp job that uses the above defined input set
        first param in VaspJob = job_cmd eg:- ["mpirun", "pvasp.5.2.11"]
        create multiple vaspjob object, each with differnt dictvaspinputset object
        consider subclassing VaspJob and overriding the post_process method
        create a custodian task using the jobs and error handlers
        consider subclassing Custodian and overriding _run_job method to modify
         when _do_check method is called
        it must be called when the calcualtion is done
        for example: check the OUTCAR and check for
        'writing wavefunctions' or maybe not

        disable error handlers on the job and the postprocess since
        the 'job' here just submitting the job to the queue 

        """
        #override the job_cmd if provided
        if job_cmd :
            for j in self.jobs:
                j.job_cmd = job_cmd
                
        c_params = {'jobs': [j.as_dict() for j in self.jobs], 'handlers': [h.as_dict() for h in self.handlers], 'max_errors': 5}
        c = Custodian(self.handlers, self.jobs, max_errors=5)
        c.run()


    def enforce_cutoff(self, input_list, delta_e=0.01):
        """
        energy difference of 10meV
        """
        matching_list = []
        for i, e in enumerate(input_list):
            if i < len(input_list)-1:
                print i, input_list[i+1], e
                print np.abs(input_list[i+1][1] - e[1])
                if np.abs(input_list[i+1][1] - e[1]) <= 0.01:
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
            print 'none of the entries satisfy the convergence criterion'
        

            
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
        drone = MPINTVaspDrone(inc_structure=True, inc_incar_n_kpoints=True) #VaspToComputedEntryDrone()#
        bg =  BorgQueen(drone)
        #bg.parallel_assimilate(rootpath)        
        bg.serial_assimilate(rootpath)
        allentries =  bg.get_data()
        alldata = []
        for e in allentries:
            if e:
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
                 setup_dir='.', parent_job_dir='.', job_dir='./Molecule'):
        
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                           setup_dir=setup_dir, parent_job_dir=parent_job_dir,
                           job_dir=job_dir)


        
    def kpoints_cnvg(self, Grid_type = 'M', kpoints_tuple = ((7, 7, 7), (10, 10, 10)),
                     conv_step = 1):
        print "Its a molecule ! no need for kpoint convergence"
        return

        
                

class CalibrateBulk(Calibrate):
    
    """
    
    Calibrate paramters for Bulk calculations
    
    """
    
    def __init__(self, incar, poscar, potcar, kpoints,
                  setup_dir='.', parent_job_dir='.', job_dir='./Bulk'):
            
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                            setup_dir=setup_dir, parent_job_dir=parent_job_dir,
                             job_dir=job_dir)

        

                    
class CalibrateSlab(Calibrate):
    
    """
    
    Calibrate paramters for Slab calculations
    
    """
    
    def __init__(self, incar, poscar, potcar, kpoints,
                 setup_dir='.', parent_job_dir='.', job_dir='./Slab'):
        
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                            setup_dir=setup_dir, parent_job_dir=parent_job_dir,
                            job_dir=job_dir)

    def kpoints_cnvg(self, Grid_type = 'M', kpoints_list = None, conv_step = 1):
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
                        print 'KPOINTmesh = ', name
                        job_dir = self.job_dir +os.sep+ 'KPOINTS' + os.sep + name
                        vis = MPINTVaspInputSet(name, self.incar, self.poscar, self.potcar,
                                                self.kpoints)
                        job = MPINTVaspJob(["pwd"], final = True, setup_dir=self.setup_dir,
                                parent_job_dir=self.parent_job_dir, job_dir=job_dir,
                                vis=vis, auto_npar=False, auto_gamma=False)
                        
                        self.jobs.append(job)
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
    calbulk = CalibrateBulk(incar, poscar, potcar, kpoints)    
    calbulk.encut_cnvg( range(400,800,100) )
    calbulk.kpoints_cnvg( kpoints_list = [ [7, 7, 7], [11, 11, 11] ] ) 
    #the job_cmd can passed to the run
    #['qsub','job_script']
    calbulk.run(['ls','-lt'])
    
    #get all data in all the directories in the provided rootfolder, here 1/
    calbulk.knob_settings('1')
    #test enforce_cutoff
    inp_list = [ ['[[2,2,4]]', 10], ['[[2,2,5]]', 9.9], ['[[2,2,6]]', 9.895], ['[[2,2,7]]', 9.888], ['[[2,2,8]]', 9.879],]
    print calbulk.enforce_cutoff(inp_list, delta_e=0.01)




#note: write the default yaml to the directory where the jobs are run. This is useful later on for comparing different job runs in that direcctory
