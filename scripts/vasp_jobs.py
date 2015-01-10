"""
This script demonstrates the usage of the module mpinterfaces/calibrate.py

Creates a sample workflow for ENCUT convergence study
"""

import os
import sys
import shutil as shu
import subprocess as sp
import socket
from math import sqrt
from collections import OrderedDict
import numpy as np

from pymatgen.matproj.rest import MPRester
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from pymatgen.core.structure import Structure
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from mpinterfaces.calibrate import CalibrateMolecule
from mpinterfaces.interface import Interface

MAPI_KEY="dwvz2XCFUEI9fJiR"

def get_struct_from_mp(formula):
    with MPRester(MAPI_KEY) as m:
        data = m.get_data(formula)
        print "\nnumber of structures matching the chemical formula "+formula+" = ", len(data)
        for d in data:
            x = {}
            x['material_id'] = str(d['material_id'])
            structure = m.get_structure_by_material_id(x['material_id'])
            return structure

#get structure from materialsproject        
strt = get_struct_from_mp('Pt')
#create slab
iface = Interface(strt, hkl=[1,1,1], min_thick=10, min_vac=10, supercell=[1,1,1])
        
#set the input files
incar_dict = {
                 'SYSTEM': 'molecule', 
                 'ENCUT': 500, 
                 'ISIF': 2, 
                 'IBRION': 2, 
                 'ALGO': 'Normal', 
                 'ISMEAR': 1, 
                 'ISPIN': 1, 
                 'EDIFF': 1e-06, 
                 'EDIFFG': -0.001, 
                 'NPAR': 8, 
                 'SIGMA': 0.1, 
                 'PREC': 'Accurate'
    }

incar = Incar.from_dict(incar_dict)
poscar = Poscar(iface, selective_dynamics = np.ones(iface.frac_coords.shape))
potcar = Potcar(poscar.site_symbols)
kpoints = Kpoints()

#set job list
encut_list = [400,500] #range(400,800,100)
turn_knobs = OrderedDict(
    [
        ('ENCUT', encut_list)
    ])

#set job command
job_dir = '/home/km468/Software/test/Molecule'
qadapter = None
job_cmd = None

nprocs = 16
incar['NPAR'] = int(sqrt(nprocs))

#on hipergator
if 'gator' in socket.gethostname():
    nnodes = 1
    walltime = '24:00:00'
    d = {'type':'PBS',
     'params':
     {
         'nnodes': str(nnodes),
         'ppnode': str(int(nprocs/nnodes)),
         'walltime': walltime,
         'job_name': 'vasp_job',
         'rocket_launch': 'mpirun ~/Software/vasp.5.3.5/vasp'
     }
    }
    qadapter = CommonAdapter(d['type'], **d['params'])
    cal = CalibrateMolecule(incar, poscar, potcar, kpoints, 
                            turn_knobs=turn_knobs, 
                            job_dir=job_dir,
                            qadapter=qadapter)

#on henniggroup machines    
else:    
    job_cmd = ['nohup',
               '/opt/openmpi_intel/bin/mpirun',
               '-n', str(nprocs),
               '/home/km468/Software/VASP/vasp.5.3.5/vasp']
    cal = CalibrateMolecule(incar, poscar, potcar, kpoints, 
                            turn_knobs=turn_knobs, 
                            job_dir=job_dir,
                            job_cmd=job_cmd, wait=False)

#setup and run        
cal.setup()
cal.run()
