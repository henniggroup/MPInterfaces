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
     database settings in the directory where it is called):

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
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket

from mpinterfaces.firetasks import MPINTCalibrateTask, MPINTMeasurementTask
from mpinterfaces import get_struct_from_mp

#---------------------------------------------------------------------
# STRUCTURES
#---------------------------------------------------------------------
structure = get_struct_from_mp('Pt', MAPI_KEY="dwvz2XCFUEI9fJiR")
#primitive --> conventional cell
sa = SpacegroupAnalyzer(structure)
structure_conventional = sa.get_conventional_standard_structure()
structure = structure_conventional.copy()
poscar = Poscar(structure)#, selective_dynamics = np.ones(iface.frac_coords.shape))
potcar = Potcar(poscar.site_symbols)

#---------------------------------------------------------------------
# INITIAL INPUTSET
#---------------------------------------------------------------------
incar_dict = { 'SYSTEM': 'Pt slab',
               'ENCUT': 500, 
                'ISIF': 2, 
                'IBRION': 2, 
                'ISMEAR': 1, 
                'EDIFF': 1e-06, 
                'NPAR': 4, 
                'SIGMA': 0.1, 
                'PREC': 'Accurate'
            }
incar = Incar.from_dict(incar_dict)
kpoints = Kpoints.automatic(20)#(80)
que  = { 'nnodes':1, 
        'nprocs':16, 
        'walltime':'24:00:00',
        'job_bin':'/home/km468/Software/VASP/vasp.5.3.5/vasp'
    }
    
turn_knobs = { 'NSW' : [100],
               'VACUUM': [10],
               'THICKNESS': [10]
             }
#---------------------------------------------------------------------
# FIRETASK 1
#
# slab 110
#---------------------------------------------------------------------
calparams1 = {}
calparams1['calibrate'] = 'CalibrateInterface'
calparams1['incar'] = incar.as_dict()
calparams1['poscar'] = poscar.as_dict()
calparams1['kpoints'] = kpoints.as_dict()
calparams1['que_params'] = que
calparams1['turn_knobs'] = turn_knobs
calparams1['system'] = {'hkl':[1,1,0], 'ligand':None}
calparams1['other_params'] = { 'job_dir':'slab_110', 'is_matrix':True}

caltask1 = MPINTCalibrateTask(calparams1)

#---------------------------------------------------------------------
# FIRETASK 2
#
#slab 111
#---------------------------------------------------------------------
calparams2 = {}
calparams2 = {k:calparams1[k] for k in calparams1.keys()}
calparams2['system'] = {'hkl':[1,1,1], 'ligand':None}
calparams2['other_params'] = {'job_dir':'slab_111', 'is_matrix':True}

caltask2 = MPINTCalibrateTask(calparams2)

#---------------------------------------------------------------------
# FIRETASK 3
#
# Measurement task
#---------------------------------------------------------------------
msrparams1 = {}
msrparams1['measurement'] = 'MeasurementInterface'
msrparams1['que_params'] = que
msrparams1['other_params'] = {'job_dir':'MSR_SLAB'}

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
if len(sys.argv)>1:
    launchpad = LaunchPad(host='localhost', port=int(sys.argv[1]), 
                          name='fireworks', username="km468", 
                          password="km468")
else:
    launchpad = LaunchPad(host='localhost', port=27017, name='fireworks',
                          username="km468", password="km468")

print('fireworks in the database before adding the workflow: \n',
      launchpad.get_fw_ids())

launchpad.add_wf(wf)

print('fireworks in the database: \n', launchpad.get_fw_ids())


##ignore, for testing purposes only
#fw1 = Firework([caltask1], name="calibrate")
#fw3 = Firework(pptask, name="post_process", parents=[fw1, fw2])
#wf = Workflow([fw1], name="mpint workflow")
#wf = Workflow([fw1, fw2, fw3], name="mpint workflow")
