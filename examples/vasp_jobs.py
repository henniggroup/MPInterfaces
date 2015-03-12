from __future__ import division, unicode_literals, print_function

"""
This script demonstrates the usage of the modules
mpinterfaces/calibrate.py and mpinterfaces/measurement.py to setup and
run vasp jobs
Note: use your own materials project key to download the required structure
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
from pymatgen.io.vaspio.vasp_input import Incar, Poscar
from pymatgen.io.vaspio.vasp_input import Potcar, Kpoints
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from mpinterfaces import get_struct_from_mp, Interface
from mpinterfaces.measurement import MeasurementInterface
from mpinterfaces.calibrate import Calibrate, CalibrateSlab
from mpinterfaces.calibrate import CalibrateInterface

#--------------------------------------------------------------------
# STRUCTURE
#--------------------------------------------------------------------
#get structure from materialsproject, use your own key       
strt = get_struct_from_mp('PbS', MAPI_KEY="dwvz2XCFUEI9fJiR")
#convert from fcc primitive to conventional cell
#the conventional unit cell is used to create the slab
#this is important becasue the hkl specification for the required slab
#is wrt the provided unit cell
sa = SpacegroupAnalyzer(strt)
structure_conventional = sa.get_conventional_standard_structure()
strt = structure_conventional.copy()
#create slab
iface = Interface(strt, hkl=[1,1,1], min_thick=10, min_vac=10,
                  supercell=[1,1,1])

#--------------------------------------------------------------------
# VASP INPUT FILES
#--------------------------------------------------------------------
incar_dict = {
                 'SYSTEM': 'Pt slab', 
                 'ENCUT': 500, 
                 'ISIF': 2, 
                 'IBRION': 2, 
                 'ISMEAR': 1, 
                 'EDIFF': 1e-06, 
                 'NPAR': 8, 
                 'SIGMA': 0.1, 
                 'NSW' : 100,
                 'PREC': 'Accurate'
    }
incar = Incar.from_dict(incar_dict)
poscar = Poscar(iface)#, selective_dynamics = np.ones(iface.frac_coords.shape))
potcar = Potcar(poscar.site_symbols)
kpoints = Kpoints.automatic(20)#(80)

#--------------------------------------------------------------------
# JOB DEFINITIONS
#--------------------------------------------------------------------
#set job list. if empty a single job will be run with the
#given inputset
encut_list = [] #range(400,800,100)
turn_knobs = OrderedDict(
    [
        ('ENCUT', encut_list)
    ])
#job directory for calibration runs
cal_job_dir = 'CAL_DIR'
#job directory for measurement runs
msr_job_dir = 'MSR_DIR'

#--------------------------------------------------------------------
# COMPUTATIONAL RESOURCE SETTINGS
#--------------------------------------------------------------------
qadapter = None
job_cmd = None
nprocs = 16
nnodes = 1
walltime = '24:00:00'
incar['NPAR'] = int(sqrt(nprocs))
d = {}
wait = True
#hipergator
if 'ufhpc' in socket.gethostname():
    d = {'type':'PBS',
     'params':
     {
         'nnodes': str(nnodes),
         'ppnode': str(int(nprocs/nnodes)),
         'walltime': walltime,
         'job_name': 'vasp_job',
         'pre_rocket': '#PBS -l pmem=1000mb',
         'rocket_launch': 'mpirun /home/km468/Software/VASP/vasp.5.3.5/vasp'
     }
    }
#stampede
elif 'stampede' in socket.gethostname():
    d = {'type':'SLURM',
     'params':
     {
         'nodes': str(nnodes),
         'ntasks': str(nprocs),
         'walltime': walltime,
         'queue':'normal',
         'account':'TG-DMR050028N',
         'job_name': 'vasp_job',
         'rocket_launch': 'ibrun /home1/01682/km468/Software/VASP/vasp.5.3.5/vasp'
     }
    }
#henniggroup machines    
elif socket.gethostname() in ['hydrogen', 'helium', 'lithium', 'beryllium', 'carbon']:    
    job_cmd = ['nohup',
               '/opt/openmpi_intel/bin/mpirun',
               '-n', str(nprocs),
               '/home/km468/Software/VASP/vasp.5.3.5/vasp']
    wait = False    
else:
    job_cmd = ['ls', '-lt']
    wait = False
    
if d:    
    qadapter = CommonAdapter(d['type'], **d['params'])

#--------------------------------------------------------------------
# setup calibration jobs and run
#--------------------------------------------------------------------
calibrations = []
cal = CalibrateSlab(incar, poscar, potcar, kpoints, system=iface.to_dict(),
                    turn_knobs=turn_knobs, qadapter=qadapter,
                    job_cmd = job_cmd, job_dir=cal_job_dir, wait=wait)
calibrations.append(cal)
for c in calibrations:
    c.setup()
    c.run()

#-------------------------------------------------------------------- 
# setup measurement jobs and run
#--------------------------------------------------------------------
measure = MeasurementInterface(calibrations, job_dir=msr_job_dir)
#CAUTION: an infinte loop.
#Breaks only if the job is done i.e the OUTCAR files from all the runs
#in the calibration jobs have the final time stamp present in it
while True:
    done  = Calibrate.check_calcs(calibrations)
    if done:
        print('calibration done ...')                
        measure.setup()
        measure.run()
        break
    else:
        print('calibrating  ...')

