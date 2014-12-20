"""
Create a workflow and submit it to the 'fireworks' database on hydrogen

Note 1:
      the database for job submission('fireworks') is
      on the mongodb server running on hydrogen.
      use own account to run the script.
      contact me(km468@cornell.edu) to create a database account
      
Note 2:
        Since hydrogen is part of the hermes subnetwork, direct connection to the database
        is not possible. So tunell port number 27017 from your local machine to port 27017 on hydrogen via ssh:
        ssh -N -f -L 27017:10.1.255.101:27017 username@hermes.mse.ufl.edu
"""

import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket #rapidfire

from mpinterfaces.firetasks import MPINTCalibrateTask, MPINTMeasurementTask

#---------------------------------------------------------
# set up the LaunchPad i.e connect to the database
#--------------------------------------------------------
launchpad = LaunchPad(host='localhost', port=27017, name='fireworks',
                       username="km468", password="km468")

#---------------------------------
# INCAR
#---------------------------------
incarparams = {'System':'test',
                   'ENCUT': 400,
                   'ISMEAR': 1,
                   'SIGMA': 0.1,
                   'EDIFF':1E-6}
incar = Incar(params=incarparams)

#---------------------------------
# POSCAR
#----------------------------------
system = 'Pt bulk'
atoms = ['Pt']
    
a0 = 3.965
lvec = [ [0.5, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5] ]
lvec = np.array(lvec) * a0
lattice = Lattice(lvec)#.from_parameters(3.866, 3.866, 3.866, 60, 60, 60)
structure = Structure( lattice, atoms, [ [0.0, 0.0, 0.0] ],coords_are_cartesian=False)
poscar = Poscar(structure, comment=system,
        selective_dynamics=None,
        true_names=True, velocities=None, predictor_corrector=None)

#---------------------------------
# KPOINTS
#----------------------------------
kpoints = Kpoints(kpts=((8, 8, 8),))

#------------------------------------------
# define firetasks
#-----------------------------------------

#-------------------------------------------------
#first firetask
#create a calibrate task: calibrate bulk
#------------------------------------------------

#set the paramters for the calibrate task
calparams1 = {}
#required params: incar, poscar, list of encut and kpoints and the calibrate object
#potcar created from poscar
calparams1['incar'] = incar.as_dict()
calparams1['poscar'] = poscar.as_dict()
calparams1['kpoints'] = kpoints.as_dict()
#submit script settings
calparams1['que'] = {
                     'type':'PBS',
                     'params':
                     {
                     'nnodes': '1',
                     'ppnode': '8',
                     'walltime': '24:00:00',
                     'job_name': 'test_job',
                     'rocket_launch': 'mpirun ~/Software/vasp.5.3.5/vasp'
                     }
                     }
#range specification for encut and kpoints
encut_list = [str(encut) for encut in range(400, 900, 100)]
kpoints_list = [str([7,7,7]) , str([11,11,11])]
#type of calibration to be done: basically the name of calibrate calss to
#be used. available options: CalibrateMolecule, CalibrateSlab, CalibrateBulk
calparams1['calibrate'] = 'CalibrateBulk'
#for other incar paramters, creat a list of values for the paramter and add it
# to the dictonary with the paramater name as the key
calparams1['turn_knobs'] = { 'ENCUT' : encut_list,
                             'KPOINTS': kpoints_list }
#optional param: job_dir is the name of the directory within which
#the calibration jobs will be run
calparams1['cal_construct_params'] = { 'job_dir':'Bulk_test'}

caltask1 = MPINTCalibrateTask(calparams1)

#---------------------------------------------------
#second firetask
#create a calibrate task: calibrate bulk calculation
#---------------------------------------------------
calparams2 = {}
calparams2 = {k:calparams1[k] for k in calparams1.keys()}
calparams2['calibrate'] = 'CalibrateSlab'
calparams2['cal_construct_params'] = {'job_dir':'Slab_1'}

caltask2 = MPINTCalibrateTask(calparams2)

#---------------------------------------------------
#third firetask
#create a Measurement task
#---------------------------------------------------
msrparams1 = {}
msrparams1['cal_objs'] = [calparams1]#, calparams2]
msrparams1['msr_construct_params'] = {'job_dir':'Measurement_1'}
msrtask1 = MPINTMeasurementTask(msrparams1)

#--------------------------------------------------
#create the fireworks from the firetasks 
#--------------------------------------------------
#fw1 = Firework([caltask1, caltask2], name="calibrate")
fw1 = Firework([caltask1], name="calibrate")
fw2 = Firework([msrtask1], name="measurement", parents=[fw1])
#fw3 = Firework(pptask, name="post_process", parents=[fw1, fw2])

#-----------------------------------------------------
#create workflow from the fireworks
#-----------------------------------------------------
#wf = Workflow([fw1], name="mpint workflow")
wf = Workflow([fw1, fw2], name="mpint workflow")
#wf = Workflow([fw1, fw2, fw3], name="mpint workflow")

print 'fireworks in the database before adding the workflow: \n', launchpad.get_fw_ids()
#-----------------------------------------
# add workflow to the database
#-----------------------------------------
launchpad.add_wf(wf)

print 'fireworks in the database: \n', launchpad.get_fw_ids()

