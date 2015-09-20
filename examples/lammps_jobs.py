from __future__ import division, unicode_literals, print_function

import os
import sys
from collections import OrderedDict

from mpinterfaces import get_struct_from_mp
from mpinterfaces.lammps import MPINTLammps, CalibrateLammps
from mpinterfaces.utils import get_run_cmmnd

# list of structures
structures = get_struct_from_mp('ZnO', all_structs=True)
# scale the structures
for s in structures:
    scell_size = 12
    a, b, c = s.lattice.abc
    s.make_supercell([int(scell_size/a),
                      int(scell_size/b),
                      int(scell_size/c)])
# lammps input paramaters    
parameters = {'atom_style': 'charge',
              'pair_style' : 'comb3 polar_off',
              'minimize':'1.0e-13  1.0e-20  1000  10000',
              'fix':['fix_nve all nve',
                     '1 all box/relax aniso 0.0 vmax 0.001',
                     '1a all qeq/comb 1 0.0001 file fq.out'] }
# pair coefficient file
pair_coeff_file = os.path.join(os.getcwd(), "ffield.comb3")
#set jobs for the list of MPINTLammps objects
turn_knobs = OrderedDict(
    [
        ('STRUCTURES', structures),
        ('PAIR_COEFFS', [pair_coeff_file])
    ])
#job directory and run settings
job_dir = 'lammps_job'
nprocs = 8
nnodes = 1
walltime = '00:30:00'
mem = 1000
job_bin = '/home/km468/Software/lammps/src/lmp_ufhpc < inp'
qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                  walltime=walltime,
                                  job_bin=job_bin, mem=mem)
# setup calibration jobs and run
cal = CalibrateLammps(parameters, turn_knobs=turn_knobs,
                      qadapter=qadapter, job_cmd = job_cmd,
                      job_dir=job_dir)
cal.setup()
cal.run()
