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
import pprint

import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar
from pymatgen.io.vaspio.vasp_input import Potcar, Kpoints
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from fireworks.fw_config import LAUNCHPAD_LOC
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
                   'EDIFF': 1e-05, 
                   'NPAR': 4, 
                   'SIGMA': 0.1, 
                   'PREC': 'Accurate'
                 }
    incar = Incar.from_dict(incar_dict)
    kpoints = Kpoints.monkhorst_automatic(kpts=(8, 8, 1))
    que  = { 'nnodes':1,
             'nprocs':16,
             'walltime':'48:00:00',
             'job_bin': '/home/km468/Software/VASP/vaspsol_kappa.5.3.5/vasp'             
            }
    # relaxation
    turn_knobs = { 'NSW' : [1000],
                   'VACUUM': [30],
                   'THICKNESS': [10]
                 }    
    is_matrix = True
    from_ase = True
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
        'is_matrix':is_matrix,
        'from_ase':from_ase,
        'Grid_type':'M'
        }
    return MPINTCalibrateTask(calparams)
    

def get_workflow(structure, hkl, wf_id=100):
    """
    returns a workflow consisting of 2 fireworks.
    firework1 has 3 calibration tasks(relax 100, 110, 111 slabs)
    firework2 has 1 measurement task(static calulations)
    """
    name = structure.composition.reduced_formula \
        +'_{0[0]}{0[1]}{0[2]}'.format(hkl)
    # calibration task1: relax hkl
    caltask = get_calibration_task(structure, hkl=hkl)
    # measurement task: static 
    msrparams = {}
    msrparams['measurement'] = 'MeasurementInterface'
    msrparams['que_params'] =  { 'nnodes':1,
                                 'nprocs':16,
                                 'walltime':'48:00:00',
                                 'job_bin': '/home/km468/Software/VASP/vaspsol_kappa.5.3.5/vasp'
                               }
    msrparams['other_params'] = {'job_dir': name+'_static'}
    msrtask = MPINTMeasurementTask(msrparams)
    # measurement task: solvation
    solmsrparams = {}
    solmsrparams['measurement'] = 'MeasurementSolvation'
    solmsrparams['que_params'] =  msrparams['que_params']
    sol_params = { 'EB_K':[78.4],
                   'TAU':[],
                   'LAMBDA_D_K':[3.0],
                   'NELECT':list(np.linspace(-1,1,5))
                 }
    solmsrparams['other_params'] = {'job_dir': name+'_sol',
                                    'sol_params':sol_params
                                   }
    solmsrtask = MPINTMeasurementTask(solmsrparams)
    #firework1 = [caltask]
    fw_calibrate = Firework([caltask],
                             name="relaxation",
                             fw_id = wf_id)
    #firework2 = [msrtask]
    msrtask['fw_id'] = wf_id+1    
    fw_measure = Firework([msrtask],
                            name="static",
                            fw_id = wf_id+1 )
    #firework3 = [solmsrtask]
    solmsrtask['fw_id'] = wf_id+2    
    fw_solmeasure = Firework([solmsrtask],
                            name="solvation",
                            fw_id = wf_id+2 )
    #workflow = [firework1, firework2, firework3]
    #firework2 is linked to firework1
    return Workflow(
        [fw_calibrate, fw_measure, fw_solmeasure],
        links_dict = {
            fw_calibrate.fw_id:[fw_measure.fw_id],
            fw_measure.fw_id:[fw_solmeasure.fw_id] } ,
        name=name )

    
if __name__=='__main__':
    species = ['Pt', 'Au', 'Ag', 'Cu']
    hkls = [[1,0,0], [1,1,0], [1,1,1]]
    workflows = []

    for i, sp in enumerate(species):
        structure = get_struct_from_mp(sp, MAPI_KEY="dwvz2XCFUEI9fJiR")
        #primitive --> conventional cell
        sa = SpacegroupAnalyzer(structure)
        structure_conventional = sa.get_conventional_standard_structure()
        structure = structure_conventional.copy()
        structure.sort()
        for j, hkl in enumerate(hkls):
            workflows.append( get_workflow(structure, hkl, wf_id=100*(i+1)+10*j) )
        
    # connect to the fireworks database and add workflow to it
    # use your own account
    launchpad = LaunchPad.from_file(LAUNCHPAD_LOC)

    print('fireworks in the database before adding the workflow: \n',
        launchpad.get_fw_ids())

    for wf in workflows:
        launchpad.add_wf(wf, reassign_all=False)

    wfs = launchpad.get_wf_ids()
    pp = pprint.PrettyPrinter(indent=4)
    for i,fw_id in enumerate(wfs):
        print('Workflow: {}'.format(i+1))
        pp.pprint(launchpad.get_wf_summary_dict(fw_id))
        print('\n')
    print('Total number of workflows in the database: \n', len(wfs))
    print('Total number of fireworks in the database: \n', len(launchpad.get_fw_ids()))
        


