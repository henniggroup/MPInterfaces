from __future__ import division, unicode_literals, print_function

"""
This is a simpler verison of single_workflow with just one firework
and firetask.
used only for testing.
Shows how to set a firework's launch directory
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

calparams1 = {}
calparams1['incar'] = incar.as_dict()
calparams1['poscar'] = poscar.as_dict()
calparams1['kpoints'] = kpoints.as_dict()
calparams1['que_params'] = { 'nnodes':1, 
                             'nprocs':16, 
                             'walltime':'24:00:00',
                             'job_bin':'/home/km468/Software/VASP/vasp.5.3.5/vasp'
                             }
turn_knobs = { 'ENCUT' : range(400, 500, 100),
               'KPOINTS': [k for k in range(20, 30, 10)]
             }
calparams1['calibrate'] = 'CalibrateBulk'
calparams1['turn_knobs'] = turn_knobs
calparams1['other_params'] = { 'job_dir':'calBulk'}

caltask1 = MPINTCalibrateTask(calparams1)

#firework with launch directory set to $FW_JOB_DIR, an environment variable
fw_calibrate = Firework([caltask1], spec={'_launch_dir':'$FW_JOB_DIR'}, name="fw_test")
wf = Workflow([fw_calibrate], name="mpint_wf_test")

lp = LaunchPad.from_file(LAUNCHPAD_LOC)

print('fireworks in the database before adding the workflow: \n',
      lp.get_fw_ids())

lp.add_wf(wf)

print('fireworks in the database: \n', lp.get_fw_ids())

