from __future__ import division, unicode_literals, print_function

"""

slab relaxation and static calculation workflows

1. get the initial structures from the materials project database
2. set relaxation calibration task for the hkl slab --> firework1
3. set static measurement task for each of the structures --> firework2
4. construct workflow consiting of the above mentioned fireworks for
   each of the structures and hkls
5. submit all the workflows to the database on hydrogen

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


def get_calibration_task(structure, hkl=[1,0,0]):
    """
    returns relaxation calibration task for hkl surface of
    the given structure
    """
    poscar = Poscar(structure)
    incar_dict = { 'SYSTEM': 'slab',
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
    # relaxation
    turn_knobs = { 'NSW' : [100],
                   'VACUUM': [5],
                   'THICKNESS': [5]
                 }    
    is_matrix = True
    # calibration task: relax hkl
    calparams = {}
    calparams['calibrate'] = 'CalibrateInterface'
    calparams['incar'] = incar.as_dict()
    calparams['poscar'] = poscar.as_dict()
    calparams['kpoints'] = kpoints.as_dict()
    calparams['que_params'] = que
    calparams['turn_knobs'] = turn_knobs
    calparams['system'] = {'hkl':hkl,
                           'ligand':None
                          }
    calparams['other_params'] = {
        'job_dir':structure.composition.reduced_formula+\
        '_{0[0]}{0[1]}{0[2]}'.format(hkl),
        'is_matrix':is_matrix
        }
    return MPINTCalibrateTask(calparams)
    

def get_workflows(structure, hkl, wf_id=100):
    """
    returns a workflow consisting of 2 fireworks.
    firework1 has 3 calibration tasks(relax 100, 110, 111 slabs)
    firework2 has 1 measurement task(static calulations)
    """
    # calibration task1: relax hkl
    caltask = get_calibration_task(structure, hkl=hkl)

    # measurement task: static 
    msrparams = {}
    msrparams['measurement'] = 'MeasurementInterface'
    msrparams['que_params'] =  { 'nnodes':1,
                                  'nprocs':16,
                                  'walltime':'24:00:00',
                                }
    msrparams['other_params'] = {
        'job_dir':structure.composition.reduced_formula \
        +'_static_measurements'
        }
    msrtask = MPINTMeasurementTask(msrparams)
    #firework1 = [caltask]
    fw_calibrate = Firework([caltask],
                             name="fw_calibrate",
                             fw_id = wf_id)
    msrtask['fw_id'] = wf_id+1
    #firework2 = [msrtask]
    fw_measure = Firework([msrtask],
                            name="fw_measurement",
                            fw_id = wf_id+1 )
    #workflow = [firework1, firework2]
    #firework2 is linked to firework1
    return Workflow( [fw_calibrate, fw_measure],
               links_dict = {fw_calibrate.fw_id:[fw_measure.fw_id]},
               name="mpint_workflow1" )

    
if __name__=='__main__':
    species = ['Pt', 'Au', 'Ag']
    hkls = [[1,0,0], [1,1,0], [1,1,1]]
    workflows = []

    for i, sp in enumerate(species):
        structure = get_struct_from_mp(sp, MAPI_KEY="dwvz2XCFUEI9fJiR")
        #primitive --> conventional cell
        sa = SpacegroupAnalyzer(structure)
        structure_conventional = sa.get_conventional_standard_structure()
        structure = structure_conventional.copy()
        structure.sort()
        for hkl in hkls:
            workflows.append( get_workflows(structure, hkl, wf_id=100*(i+1)) )
        
    # connect to the fireworks database and add workflow to it
    # use your own account
    if len(sys.argv)>1:
        launchpad = LaunchPad(host='localhost', port=int(sys.argv[1]),
                            name='fireworks', username="km468", 
                            password="km468" )
    else:
        launchpad = LaunchPad(host='localhost', port=27017,
                            name='fireworks', username="km468",
                            password="km468" )

    print('fireworks in the database before adding the workflow: \n',
        launchpad.get_fw_ids())

    for wf in workflows:
        launchpad.add_wf(wf, reassign_all=False)

    print('fireworks in the database: \n', launchpad.get_fw_ids())


