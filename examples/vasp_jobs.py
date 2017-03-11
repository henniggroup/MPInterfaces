# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This script demonstrates the usage of the module
mpinterfaces/calibrate.py to setup and run vasp jobs
"""

import os
from math import sqrt
from collections import OrderedDict

from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.inputs import Potcar, Kpoints
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from mpinterfaces import get_struct_from_mp
from mpinterfaces.calibrate import CalibrateSlab
from mpinterfaces.interface import Interface
from mpinterfaces.utils import get_run_cmmnd

MAPI_KEY = os.environ.get("MAPI_KEY", "")
# get structure from materialsproject, use your own key
strt = get_struct_from_mp('PbS', MAPI_KEY=MAPI_KEY)
# convert from fcc primitive to conventional cell
# the conventional unit cell is used to create the slab
# this is important becasue the hkl specification for the required slab
# is wrt the provided unit cell
sa = SpacegroupAnalyzer(strt)
structure_conventional = sa.get_conventional_standard_structure()
strt = structure_conventional.copy()
# create slab
iface = Interface(strt, hkl=[1, 1, 1],
                  min_thick=10, min_vac=10,
                  supercell=[1, 1, 1])
# sort structure into groups of elements atoms for Potcar mapping 
iface.sort()
# vasp input
incar_dict = {
    'SYSTEM': 'test',
    'ENCUT': 500,
    'ISIF': 2,
    'IBRION': 2,
    'ISMEAR': 1,
    'EDIFF': 1e-06,
    'NPAR': 8,
    'SIGMA': 0.1,
    'NSW': 100,
    'PREC': 'Accurate'
}
incar = Incar.from_dict(incar_dict)
poscar = Poscar(iface)
potcar = Potcar(poscar.site_symbols)
kpoints = Kpoints.automatic(20)
# set job list. if empty a single job will be run
encut_list = []  # range(400,800,100)
turn_knobs = OrderedDict(
    [
        ('ENCUT', encut_list)
    ])
# job directory
job_dir = 'vasp_job'
# run settings
qadapter = None
job_cmd = None
nprocs = 16
nnodes = 1
walltime = '24:00:00'
mem = 1000
incar['NPAR'] = int(sqrt(nprocs))
job_bin = '/home/km468/Software/VASP/vasp.5.3.5/vasp'
qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                  walltime=walltime,
                                  job_bin=job_bin, mem=mem)
# setup calibration jobs and run
cal = CalibrateSlab(incar, poscar, potcar, kpoints,
                    system=iface.as_dict(),
                    turn_knobs=turn_knobs, qadapter=qadapter,
                    job_cmd=job_cmd, job_dir=job_dir)
cal.setup()
cal.run()
