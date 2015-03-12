from __future__ import division, unicode_literals, print_function

"""
Create a solvation measurement firework from the cal_objs of an
already existing firework in the fireworks database

CAUTION: Do not run the script as such
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


def get_workflow(name, lp, fw_id):
    """
    measure firework from the cal_objs of firework with id fw_id
    """
    fw = lp.get_fw_by_id(fw_id)
    que_params =  { 'nnodes':1,
                    'nprocs':16,
                    'walltime':'24:00:00',
                    'job_bin': '/home/km468/Software/VASP/vaspsol_kappa.5.3.5/vasp'
                  }
    # measurement task: solvation
    solmsrparams = {}
    solmsrparams['measurement'] = 'MeasurementSolvation'
    solmsrparams['que_params'] =  que_params
    sol_params = { 'EB_K':[78.4],
                   'TAU':[],
                   'LAMBDA_D_K':[3.0],
                   'NELECT':[0]
                 }
    solmsrparams['other_params'] = {'job_dir': name,
                                    'sol_params':sol_params
                                   }
    solmsrtask = MPINTMeasurementTask(solmsrparams)
    ##firework = [solmsrtask]
    spec = {}
    spec['cal_objs'] = fw.spec['cal_objs']
    #spec['_launch_dir'] = '/scratch/lfs/km468/MP/kappa/30'
    #spec['_launch_dir'] = '$FW_JOB_DIR/KAPPA'
    fw_solmeasure = Firework([solmsrtask], spec=spec, name="solvation_pbz")
    return Workflow([fw_solmeasure],name=name)

    
if __name__=='__main__':
    lp = LaunchPad.from_file(LAUNCHPAD_LOC)

    print('fireworks in the database(before): \n', lp.get_fw_ids())
            
    species = ['Pt', 'Au', 'Ag', 'Cu'] #100, 200, 300, 400
    hkls = [[1,0,0], [1,1,0], [1,1,1]]
    workflows = []

    for i,sp in enumerate(species):
        for j, hkl in enumerate(hkls):
            fw_id = 666
            name = sp+'_{0[0]}{0[1]}{0[2]}'.format(hkl)+'_npbz'
            workflows.append(get_workflow(name, lp, fw_id))
            
    for wf in workflows:
        lp.add_wf(wf, reassign_all=False)
