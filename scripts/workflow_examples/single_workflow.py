from __future__ import division, unicode_literals, print_function

"""
Create a workflow and submit it to the 'fireworks' database on hydrogen

Note 1:
      the database for job submission('fireworks') is
      on the mongodb server running on hydrogen.
      use own account to run the script.
      contact me(km468@cornell.edu) to create a database account
      
Note 2:
        Since hydrogen is part of the hermes subnetwork, direct
        connection to the database is not possible. So tunell 
        port number 27017 from your local machine to port 27017 
        on hydrogen via ssh:

        ssh -N -f -L 27017:10.1.255.101:27017 username@hermes.mse.ufl.edu

        if port 27017 on the machine that you are running is 
        not available, use another port number for tunneling. example:-

        ssh -N -f -L 27030:10.1.255.101:27017 username@hermes.mse.ufl.edu

        mind: if the tunneled port is changed, the port number
        in the launchpad initialization should also be changed

Note 3:
     to submit workflow to the database:

         python submit_workflow.py
 
     Fireworks package has some nice utility scripts for 
     launching fireworks and checking job status. If the 
     fireworks package is installed then those scripts are already
     in your PATH. Some examples are given below

     initialize database connection(writes a yaml file with the 
     database settings in the directory where it is called).
     Tip: dont have to do this again if you create a ~/.fireworks
     folder and the generated yaml file there:

         lpad init

     CAUTION: Be careful when using the following command as it will 
     erase all workflows from the database:

         lpad reset

     launch a single firework:

         rlaunch singleshot

     get all fireworks info:

         lpad get_fws
"""

import sys

import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints

from fireworks.fw_config import LAUNCHPAD_LOC
from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket

from mpinterfaces.firetasks import MPINTCalibrateTask, MPINTMeasurementTask

#---------------------------------------------------------------------
# INITIAL INPUTSET
#---------------------------------------------------------------------
#structure
a0 = 3.965
lattice_matrix = np.array([ [0.5, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5] ]) * a0
lattice = Lattice(lattice_matrix)
structure = Structure( lattice, ['Pt'], [ [0.0, 0.0, 0.0] ],
                       coords_are_cartesian=False)
incarparams = {'System':'test',
                   'ENCUT': 400,
                   'ISMEAR': 1,
                   'SIGMA': 0.1,
                   'EDIFF':1E-6}
incar = Incar(params=incarparams)
poscar = Poscar(structure, comment='test', selective_dynamics=None)
potcar = Potcar(symbols = poscar.site_symbols, functional='PBE',
                sym_potcar_map=None)
kpoints = Kpoints(kpts=((8, 8, 8),))

#---------------------------------------------------------------------
# FIRETASK
#
# calibratebulk task
#---------------------------------------------------------------------
calparams1 = {}
calparams1['incar'] = incar.as_dict()
calparams1['poscar'] = poscar.as_dict()
calparams1['kpoints'] = kpoints.as_dict()
calparams1['que_params'] = { 'nnodes':1, 
                             'nprocs':16, 
                             'walltime':'24:00:00',
                             'job_bin':'/home/km468/Software/VASP/vasp.5.3.5/vasp'
                             }
#calparams1['job_cmd'] = job_cmd
turn_knobs = { 'ENCUT' : range(400, 900, 100),
               'KPOINTS': [k for k in range(20, 40, 10)]
             }
#type of calibration to be done: basically the name of calibrate class
#to be used. 
calparams1['calibrate'] = 'CalibrateBulk'
calparams1['turn_knobs'] = turn_knobs
calparams1['other_params'] = { 'job_dir':'calBulk'}

caltask1 = MPINTCalibrateTask(calparams1)

#---------------------------------------------------------------------
# FIRETASK
#
# calibrateinterface task
#---------------------------------------------------------------------
calparams2 = {}
calparams2 = {k:calparams1[k] for k in calparams1.keys()}
calparams2['calibrate'] = 'CalibrateInterface'
calparams2['system'] = {'hkl':[1,1,1], 'ligand':None}
calparams2['other_params'] = {'job_dir':'califace'}

caltask2 = MPINTCalibrateTask(calparams2)

#---------------------------------------------------------------------
# FIRETASK
#
# Measurement task
#---------------------------------------------------------------------
msrparams1 = {}
msrparams1['measurement'] = 'MeasurementInterface'
msrparams1['que_params'] = calparams1['que_params']
msrparams1['other_params'] = {'job_dir':'Measurement_1'}

msrtask1 = MPINTMeasurementTask(msrparams1)

#---------------------------------------------------------------------
# FIREWORKS
#
# create fireworks from the firetasks
#---------------------------------------------------------------------
fw_calibrate = Firework([caltask1, caltask2], name="fw_calibrate")
fw_measure = Firework([msrtask1], name="fw_measurement", parents=[fw_calibrate])

#---------------------------------------------------------------------
# WORKFLOW
#
# create workflow from the fireworks
#---------------------------------------------------------------------
wf = Workflow([fw_calibrate, fw_measure], name="mpint_workflow")

#---------------------------------------------------------------------
# connect to the fireworks database and add workflow to it
# use your own account
#---------------------------------------------------------------------
launchpad = LaunchPad.from_file(LAUNCHPAD_LOC)

print('fireworks in the database before adding the workflow: \n',
      launchpad.get_fw_ids())

launchpad.add_wf(wf)

print('fireworks in the database: \n', launchpad.get_fw_ids())


##ignore, for testing purposes only
#fw1 = Firework([caltask1], name="calibrate")
#fw3 = Firework(pptask, name="post_process", parents=[fw1, fw2])
#wf = Workflow([fw1], name="mpint workflow")
#wf = Workflow([fw1, fw2, fw3], name="mpint workflow")
