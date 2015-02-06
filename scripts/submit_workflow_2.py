from __future__ import division, unicode_literals, print_function

"""
Create and submit multiple(3) workflows it to the 'fireworks' database
on hydrogen

"""

import sys

import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar
from pymatgen.io.vaspio.vasp_input import Potcar, Kpoints
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket

from mpinterfaces import get_struct_from_mp
from mpinterfaces.firetasks import MPINTCalibrateTask
from mpinterfaces.firetasks import MPINTMeasurementTask
from mpinterfaces.firetasks import MPINTDatabaseTask

#---------------------------------------------------------------------
# STRUCTURES
#---------------------------------------------------------------------
structure = get_struct_from_mp('Pt', MAPI_KEY="dwvz2XCFUEI9fJiR")
#primitive --> conventional cell
sa = SpacegroupAnalyzer(structure)
structure_conventional = sa.get_conventional_standard_structure()
structure = structure_conventional.copy()
poscar = Poscar(structure)
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
       }
turn_knobs = { 'NSW' : [100],
               'VACUUM': [5],
               'THICKNESS': [5]
             }    
is_matrix = True

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
calparams1['system'] = {'hkl':[1,1,0],
                        'ligand':None
                       }
calparams1['other_params'] = { 'job_dir':'slab_110',
                               'is_matrix':is_matrix
                             }

caltask1 = MPINTCalibrateTask(calparams1)

#---------------------------------------------------------------------
# FIRETASK 2
#
#slab 111
#---------------------------------------------------------------------
calparams2 = {k:calparams1[k] for k in calparams1.keys()}
calparams2['system'] = {'hkl':[1,1,1],
                        'ligand':None
                       }
calparams2['other_params'] = {'job_dir':'slab_111',
                              'is_matrix':is_matrix
                             }

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
fw_calibrate1 = Firework( [caltask1, caltask2],
                         name="fw_calibrate",
                         fw_id = 10 )
fw_measure1 = Firework( [msrtask1],
                       name="fw_measurement",
                       fw_id = 11 )

fw_calibrate2 = Firework( [caltask1, caltask2],
                         name="fw_calibrate",
                         fw_id = 20 )
fw_measure2 = Firework( [msrtask1],
                       name="fw_measurement",
                       fw_id = 21 )

fw_calibrate3 = Firework( [caltask1, caltask2],
                         name="fw_calibrate",
                         fw_id = 30 )
fw_measure3 = Firework( [msrtask1],
                       name="fw_measurement",
                       fw_id = 31 )


#---------------------------------------------------------------------
# WORKFLOW
#
# create workflow from the fireworks
#---------------------------------------------------------------------
wf1 = Workflow( [fw_calibrate1, fw_measure1],
               links_dict = {10:[11]},
               name="mpint_workflow1" )

wf2 = Workflow( [fw_calibrate2, fw_measure2],
               links_dict = {20:[21]},
               name="mpint_workflow2" )

wf3 = Workflow( [fw_calibrate3, fw_measure3],
               links_dict = {30:[31]},
               name="mpint_workflow3" )

workflows = [wf1, wf2, wf3]

#---------------------------------------------------------------------
# connect to the fireworks database and add workflow to it
# use your own account
#---------------------------------------------------------------------
if len(sys.argv)>1:
    launchpad = LaunchPad( host='localhost', port=int(sys.argv[1]), 
                           name='fireworks', username="km468", 
                           password="km468" )
else:
    launchpad = LaunchPad( host='localhost', port=27017,
                           name='fireworks', username="km468",
                           password="km468" )

print('fireworks in the database before adding the workflow: \n',
      launchpad.get_fw_ids())

for wf in workflows:
    launchpad.add_wf(wf, reassign_all=False)

print('fireworks in the database: \n', launchpad.get_fw_ids())
