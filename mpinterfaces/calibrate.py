"""
Calibration module

"""
import sys
import os
import shutil
import shlex
import subprocess
import operator
from collections import Counter, OrderedDict
import re
import time
import datetime
from pprint import pprint
import logging

import numpy as np
from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar
from pymatgen.io.vaspio.vasp_input import Potcar, Kpoints
#from custodian.vasp.handlers import VaspErrorHandler
#from custodian.vasp.handlers import FrozenJobErrorHandler
#from custodian.vasp.handlers import MeshSymmetryErrorHandler
#from custodian.vasp.handlers import NonConvergingErrorHandler
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
                 qadapter=None, job_cmd='qsub',
                 turn_knobs={'ENCUT':[],'KPOINTS':[]} ):
        """
        setup_dir = directory from where the setup files are
        copied from.
        parent_job_dir = the directory from which the script is run
        job_dir = name of the job directory relative to
        the parent directory
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
        #example:- handlers = [VaspErrorHandler(),
        #FrozenJobErrorHandler(),
        #MeshSymmetryErrorHandler(), NonConvergingErrorHandler()]
        self.handlers = [ ]
        self.job_cmd = job_cmd
        self.n_atoms = 0
        self.turn_knobs = turn_knobs
        self.response_to_knobs = {}
        self.sorted_response_to_knobs = {}        
        for k, v in turn_knobs.items():
            self.response_to_knobs[k] = {}
            self.sorted_response_to_knobs[k] = {}            

    def setup(self, conv_step = 1, matrix = False):
        if matrix == True: 
		for a, b in self.turn_knobs.items():
			if a == 'KPOINTS' and b:
				for c, d in self.turn_knobs.items():
					if c == 'ENCUT' or c == 'SIGMA':
						for e, f in self.turn_knobs.items():
							if e == 'VOLUME' and f:
								self.setup_poscar_jobs(scale_list = f)
						self.setup_incar_jobs(c, d)
		 
		
		
	else:
		for k, v in self.turn_knobs.items():
            		if k == 'KPOINTS' and v:
	                	self.setup_kpoints_jobs(kpoints_list = v, conv_step = conv_step)
            		if k == 'ENCUT' or k == 'SIGMA':
                		self.setup_incar_jobs(k, v)
	    		if k == 'VOLUME' and v:
				self.setup_poscar_jobs(scale_list = v)
            
    def setup_matrix_job(job_1_list, job_2_list):
        #for doing  matrix convergence job implementing for 2 param matrix first
        #for k, v in self.turn_knobs
	pass

    def setup_incar_jobs(self, param, val_list):
        for val in val_list:
            print 'setting INCAR parameter '+param+' = ', val
            job_dir  = self.job_dir+ os.sep + \
              param + os.sep + str(val)
            self.incar[param] = val
            self.add_job(name=param+str(val), job_dir=job_dir)
            
    def setup_kpoints_jobs(self, Grid_type = 'A',
                           kpoints_list = None, conv_step = 1):
        """
        set the jobs for kpoint convergence
        
        """
        if Grid_type == 'M':
	#Monkhorst_pack method for bulk
            	if kpoints_list:
                	conv_list = kpoints_list 
                	start = conv_list[0]
                	end = conv_list[1]
                	if (conv_step):
                    		for x in range(conv_step+start[0], end[0], conv_step):
                        		conv_list.append([x, x, x])
                    		for kpoint in conv_list:
                        		self.kpoints = \
                          		Kpoints.monkhorst_automatic(kpts = kpoint)
                        		name = str(kpoint[0]) + 'x' + \
                          		str(kpoint[1]) + 'x' + str(kpoint[2])
                        		print 'KPOINTS = ', name
                        		job_dir = self.job_dir +os.sep+ 'KPOINTS' +\
                           		os.sep + name
                        		self.add_job(name=name, job_dir=job_dir) 
            	else:
                	print 'kpoints_list not provided'
	if Grid_type == 'A':
	#Autoomatic method default 
		 if kpoints_list:
            		conv_list = kpoints_list
            		start = conv_list[0]
            		end = conv_list[1]
            	 	if (conv_step):
                		for x in range(conv_step+start, end, conv_step):
                    			conv_list.append(x)
                		for kpoint in conv_list:
                    			self.kpoints = Kpoints.automatic(subdivisions = kpoint)
                    			name = str(kpoint)
                    			print 'KPOINTS = ', name
                    			job_dir = self.job_dir +os.sep+ 'KPOINTS' +\
                    			os.sep + name
                    			self.add_job(name=name, job_dir=job_dir)
	    	 else:
        		print 'kpoints_list not provided'

		
	   
    def setup_poscar_jobs(self, scale_list, Name = "volume_scale_") :
	#for scaling the latice vectors of the original structure, scale_list is volume scaling factor list 
	print "setting volume as given list "
	Store_struct = self.poscar.structure
	Store_struct.to(fmt = "poscar", filename= "POSCAR.orig")
	for s in scale_list:
    		Scale_struct = Structure.from_file("POSCAR.orig")
    		Scale_struct.scale_lattice(s*Scale_struct.volume)
    		#print Scale_struct
    		Scale_struct.to(fmt="poscar", filename= "POSCAR")
    		self.poscar = Poscar(Scale_struct)
    		#print New_poscar
		job_dir  = self.job_dir+ os.sep + 'VOLUME' +\
              	  os.sep + str(s)
            	self.add_job(name=Name+str(s), job_dir=job_dir)

    def add_job(self, name='noname', job_dir='.'):
        vis = MPINTVaspInputSet(name, self.incar, self.poscar,
                                self.potcar, self.kpoints,
                                self.qadapter)
        #the job command can be overrridden in the run method
        job = MPINTVaspJob(self.job_cmd, final = True,
                           setup_dir=self.setup_dir,
                           parent_job_dir=self.parent_job_dir,
                           job_dir=job_dir, vis=vis, auto_npar=False,
                            auto_gamma=False)
        self.jobs.append(job)

    
    def run(self, job_cmd=None):
        """
        run the vasp jobs through custodian

        """
        #if the job list is empty,
        #run a single job with the provided input set
        if not self.jobs :
            self.add_job(name='single job', job_dir=self.job_dir)  
        
        #override the job_cmd if provided
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
        
    def set_knob_responses(self):
        drone = \
          MPINTVaspDrone(inc_structure=True, inc_incar_n_kpoints=True)
        bg =  BorgQueen(drone)
        for k, v in self.response_to_knobs.items():
            rootpath = self.job_dir+ os.sep + k
            print 'rootpath = ', rootpath
            #bg.parallel_assimilate(rootpath)        
            bg.serial_assimilate(rootpath)
            allentries =  bg.get_data()
            for e in allentries:
                if e:
                    self.n_atoms = len(e.structure)
                    #print 'n_atoms', self.n_atoms
                    if k == 'KPOINTS':
                        self.response_to_knobs[k][str(e.kpoints.kpts)] \
                           = e.energy
                    else:
                        self.response_to_knobs[k][str(e.incar[k])] \
                          = e.energy

    def set_sorted_optimum_params(self):
        matching_knob_responses = []
        sorted_knob_responses = []                
        self.optimum_knob_responses = {}
        #order the keys(encut or kpoint)from large value of the
        # energy  to small value
        for k, v in self.response_to_knobs.items():
            sorted_knob_responses = \
              sorted(v.items(), key=operator.itemgetter(1),
                     reverse=True)
            #print 'sorted_knob_responses ', sorted_knob_responses
            #get the list of encut and kpoints that
            #satisfy the delate criterion
            #mind: default deltae = 0.001eV per atom
            matching_knob_responses = \
              self.enforce_cutoff(sorted_knob_responses)
            self.sorted_response_to_knobs[k] = \
              OrderedDict(sorted_knob_responses)
            #print 'matching_knob_response ', matching_knob_responses
            if matching_knob_responses:
                if k == "KPOINTS":
                    nkpt = matching_knob_responses[0][0] * \
                      matching_knob_responses[0][1] * matching_knob_responses[0][2]
                    for i, val in enumerate(matching_knob_responses):
                        if i < len(matching_kpt)-1:
                            nkpt1 = val[i+1][0] * val[i+1][1] * \
                              val[i+1][2]
                            if nkpt1<nktp:
                                self.optimum_knob_responses[k] = val
                else:
                    self.optimum_knob_responses[k] = \
                      min(matching_knob_responses)
                        
    def enforce_cutoff(self, input_list, delta_e_peratom=0.001):
        """
        energy difference of 1meV per atom
        """
        matching_list = []
        for i, e in enumerate(input_list):
            if i < len(input_list)-1:
                #print i, input_list[i+1], e
                #print np.abs(input_list[i+1][1] - e[1])
                if np.abs(input_list[i+1][1] - e[1])/self.n_atoms <= \
                  delta_e_peratom:
                    matching_list.append(input_list[i+1][0])
        if matching_list:
            #print matching_list
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
                return [float(encut) for encut in matching_list]
        else:
            return []
                    
        
class CalibrateMolecule(Calibrate):
    """
    
    Calibrate paramters for Molecule calculations
    
    """
    def __init__(self, incar, poscar, potcar, kpoints,
                 setup_dir='.', parent_job_dir='.',
                 job_dir='./Molecule', qadapter=None,
                 job_cmd='qsub',
                 turn_knobs={'ENCUT':[],'KPOINTS':[]}):
        
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                           setup_dir=setup_dir,
                           parent_job_dir=parent_job_dir,
                           job_dir=job_dir, qadapter=qadapter,
                           job_cmd=job_cmd,
                           turn_knobs = turn_knobs)
        
    def setup_kpoints_jobs(self, Grid_type = 'M',
                           kpoints_list = None, conv_step = 1):
        print "Its a molecule ! no need for kpoint convergence"
        self.kpoints = Kpoints.monkhorst_automatic(kpts = [1,1,1])
        return

    
class CalibrateBulk(Calibrate):
    """
    
    Calibrate paramters for Bulk calculations
    
    """
    def __init__(self, incar, poscar, potcar, kpoints,
                  setup_dir='.', parent_job_dir='.',
                  job_dir='./Bulk', qadapter=None,
                  job_cmd='qsub',
                  turn_knobs={'ENCUT':[],'KPOINTS':[]}): 
            
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                            setup_dir=setup_dir,
                             parent_job_dir=parent_job_dir,
                             job_dir=job_dir, qadapter=qadapter,
                              job_cmd=job_cmd,
                             turn_knobs = turn_knobs )
        

class CalibrateSlab(Calibrate):
    """
    
    Calibrate paramters for Slab calculations
    
    """
    def __init__(self, incar, poscar, potcar, kpoints,
                 setup_dir='.', parent_job_dir='.', job_dir='./Slab',
                qadapter=None, job_cmd='qsub',
                turn_knobs={'ENCUT':[],'KPOINTS':[]}):
        
        Calibrate.__init__(self, incar, poscar, potcar, kpoints,
                            setup_dir=setup_dir,
                             parent_job_dir=parent_job_dir,
                            job_dir=job_dir, qadapter=qadapter,
                             job_cmd=job_cmd,
                            turn_knobs = turn_knobs)

    def setup_kpoints_jobs(self, Grid_type = 'M',
                            kpoints_list = None, conv_step = 1):
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
                        self.kpoints = \
                          Kpoints.monkhorst_automatic(kpts = kpoint)
                        name = str(kpoint[0]) + 'x' + \
                          str(kpoint[1]) + 'x' + str(kpoint[2])
                        print 'KPOINTS = ', name
                        job_dir = self.job_dir +os.sep+ 'KPOINTS' + \
                          os.sep + name
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
    atoms = poscar.site_symbols
    potcar = Potcar(symbols=atoms, functional='PBE',
                    sym_potcar_map=None)
    kpoints = Kpoints.monkhorst_automatic(kpts=(16, 16, 16),
                                          shift=(0, 0, 0))

    calbulk = CalibrateBulk(incar, poscar, potcar, kpoints,
                            job_dir='./Bulk',
                            turn_knobs = {'ENCUT':range(400,800,100),
                                          'KPOINTS':[
                                              20, 60
                                              ], 'VOLUME':[0.6, 0.8, 1.0, 1.2, 1.4] } )    
    calbulk.setup(conv_step = 1)
    calbulk.run(['ls','-lt'])
    #get the knob responses
    #calbulk.set_knob_responses()
    #print calbulk.response_to_knobs
    #optimu knob responses
    #calbulk.set_sorted_optimum_params()
    #print calbulk.sorted_response_to_knobs['ENCUT']['600.0']    
    #print calbulk.optimum_knob_responses
    #test enforce_cutoff
    #inp_list = [ ['[[2,2,4]]', 10], ['[[2,2,5]]', 9.9],
    #           ['[[2,2,6]]', 9.895], ['[[2,2,7]]', 9.888],
    #['[[2,2,8]]', 9.879],]
    #print calbulk.enforce_cutoff(inp_list, delta_e=0.01)




#note: write the default yaml to the directory where the jobs are run.
# This is useful later on for comparing different job runs in that
#direcctory
