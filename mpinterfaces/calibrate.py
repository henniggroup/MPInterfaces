"""
Calibration module

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
#from custodian.vasp.handlers import VaspErrorHandler, FrozenJobErrorHandler, MeshSymmetryErrorHandler, NonConvergingErrorHandler
from custodian.custodian import Custodian, gzip_dir
from custodian.vasp.interpreter import VaspModder
from mpinterfaces.instrument import MPINTVaspInputSet, MPINTVaspJob


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


            
    def kpoints_cnvg(self, Grid_type = 'M', kpoints_tuple = ((7, 7, 7), (10, 10, 10)),
                     conv_step = 1):
        """
        must be implemented in the deriving class
        """
        pass


    
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


        
    def knob_settings(self):
       """
       go through the run directories and get the optimal values for the input paramters
       """
       pass


        
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
        """
        If molecule constructs single kpoint file with MP 1x1x1.
        Calls Kpoints constructors according to Grid_type: Monkhorst Pack Automatic is default,
        G for Gamma centered automatic; ONLY MP method implemented now others can be added 
        if_slab flag is for the application of constraint to z-axis for slab,
        kpoints_list describes the start and end of kpoints set
        eg: user can pass ((6, 6, 6), ((10, 10, 10))) and conv_step for a
        convergence to be done for 6x6x6 to 10x10x10,
        if_slab can be switched to True constraints z to default to 1,
        that is 6x6x1 to 10x10x1
        first defines the list of kpoints according to the user input,
        user needs to give only the start and end kpoint,
        whether it is for a slab, molecule or simple bulk
        """
        kpoint = (1, 1, 1)
        self.kpoints = self.kpoints.monkhorst_automatic(kpts=kpoint)
        K = list(kpoint)
        print 'KPOINTmesh = ', str(K[0])+'x'+str(K[1])+'x'+str(K[2])
        job_dir = self.job_dir + os.sep + 'KPOINTS' + os.sep + str(K[0])+'x'+str(K[1])+'x'+str(K[2])
        vis = MPINTVaspInputSet('kpoint_'+str(K[0])+'x'+str(K[1])+'x'+str(K[2]),
                                self.incar, self.poscar, self.potcar, self.kpoints)
        job = MPINTVaspJob(["pwd"], final = True, setup_dir=self.setup_dir,
                               parent_job_dir=self.parent_job_dir, job_dir=job_dir,
                               vis=vis, auto_npar=False, auto_gamma=False)
        
        job = MPINTVaspJob(["pwd"], final = True,
                           setup_dir=self.setup_dir,
                           job_dir=os.path.abspath(job_dir), vis=vis,
                           auto_npar=False, auto_gamma=False)
        self.jobs.append(job)
        
                

class CalibrateBulk(Calibrate):
    
    """
    
    Calibrate paramters for Bulk calculations
    
    """
    
    def __init__(self, incar, poscar, potcar, kpoints,
                  setup_dir='.', parent_job_dir='.', job_dir='./Bulk'):
            
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                            setup_dir=setup_dir, parent_job_dir=parent_job_dir,
                             job_dir=job_dir)

        
    def kpoints_cnvg(self, Grid_type = 'M', kpoints_tuple = ((7, 7, 7), (10, 10, 10)),
                      conv_step = 1):
        if Grid_type == 'M':
            print "at Grid"
            #local list convergence_list , convert from tuple
            #because constructor takes tuple as argument            
            conv_list = list(kpoints_tuple) 
            start = list(conv_list[0])
            end = list(conv_list[1])
            if (conv_step):
                for x in range(1+start[0], end[0], conv_step):
                    conv_list.append([x, x, x])                
                for kpoint in conv_list:
                    self.kpoints = self.kpoints.monkhorst_automatic(kpts = kpoint)
                    K = list(kpoint)
                    print 'KPOINTmesh = ', str(K[0])+'x'+str(K[1])+'x'+str(K[2])
                    job_dir = self.job_dir +os.sep+ 'KPOINTS' + os.sep + str(K[0])+'x'+str(K[1])+'x'+str(K[2])
                    vis = MPINTVaspInputSet('kpoint_'+str(K[0])+'x'+str(K[1])+'x'+str(K[2]),
                                            self.incar, self.poscar, self.potcar,
                                            self.kpoints)
                    job = MPINTVaspJob(["pwd"], final = True, setup_dir=self.setup_dir,
                               parent_job_dir=self.parent_job_dir, job_dir=job_dir,
                               vis=vis, auto_npar=False, auto_gamma=False)
                    
                    self.jobs.append(job)
        

                    
class CalibrateSlab(Calibrate):
    
    """
    
    Calibrate paramters for Slab calculations
    
    """
    
    def __init__(self, incar, poscar, potcar, kpoints,
                 setup_dir='.', parent_job_dir='.', job_dir='./Slab'):
        
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                            setup_dir=setup_dir, parent_job_dir=parent_job_dir,
                            job_dir=job_dir)

        
    def kpoints_cnvg(self, Grid_type = 'M', kpoints_tuple = ((7, 7, 7), (10, 10, 10)),
                      conv_step = 1):
        if Grid_type == 'M':
            conv_list = list(kpoints_tuple) 
            start = list(conv_list[0])
            end = list(conv_list[1])
            if start[2] != 1 or end[2] != 1:
                print "Kpoints not for slab input!"
            elif (conv_step):
                for x in range(1+start[0], end[0], conv_step):
                    conv_list.append([x, x, 1])
                for kpoint in conv_list:
                    self.kpoints = self.kpoints.monkhorst_automatic(kpts = kpoint)
                    K = list(kpoint)
                    print 'KPOINTmesh = ', str(K[0])+'x'+str(K[1])+'x'+str(K[2])
                    job_dir = self.job_dir +os.sep+ 'KPOINTS' + os.sep + str(K[0])+'x'+str(K[1])+'x'+str(K[2])
                    vis = MPINTVaspInputSet('kpoint_'+str(K[0])+'x'+str(K[1])+'x'+str(K[2]),
                                            self.incar, self.poscar, self.potcar,
                                            self.kpoints)
                    job = MPINTVaspJob(["pwd"], final = True, setup_dir=self.setup_dir,
                               parent_job_dir=self.parent_job_dir, job_dir=job_dir,
                               vis=vis, auto_npar=False, auto_gamma=False)

                    self.jobs.append(job)
	            


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
    calbulk.encut_cnvg(range(400,800,100))
    calbulk.kpoints_cnvg(kpoints_tuple = ((7, 7, 7), (11, 11, 11))) 
    #the job_cmd can passed to the run
    #['qsub','job_script']
    calbulk.run(['ls','-lt'])



#note: write the default yaml to the directory where the jobs are run. This is useful later on for comparing different job runs in that direcctory
