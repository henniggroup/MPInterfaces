from __future__ import division, unicode_literals, print_function

import os
import sys
from collections import OrderedDict

from mpinterfaces import get_struct_from_mp
from mpinterfaces.lammps import CalibrateLammps
from mpinterfaces.utils import *

# all the info/warnings/outputs redirected to the log file: lammps.log
logger = get_logger('lammps')

# list of structures from materialsproject
structures = get_struct_from_mp('ZnO', all_structs=True)
# scale the structures
scell_size = 12
for s in structures:
    a, b, c = s.lattice.abc
    s.make_supercell([int(scell_size/a),
                      int(scell_size/b),
                      int(scell_size/c)])
# lammps input paramaters    
parameters = {'atom_style': 'charge',
              'charges': {'Zn':2, 'O':-2},
              'minimize':'1.0e-13  1.0e-20  1000  10000',
              'fix':['fix_nve all nve',
                     '1 all box/relax aniso 0.0 vmax 0.001',
                     '1a all qeq/comb 1 0.0001 file fq.out'] }
# list of pair styles
pair_styles = ['comb3 polar_off']
# list of pair coefficient files
pair_coeff_files = [os.path.join(os.getcwd(), "ffield.comb3")]

def test(**kwargs):
    turn_knobs = OrderedDict(
        [
            ('STRUCTURES', structures),
            ('PAIR_STYLES', pair_styles),
            ('PAIR_COEFFS', pair_coeff_files)
            ] )
    # job directory and run settings
    job_dir = 'lammps_job'
    nprocs = 4
    nnodes = 1
    walltime = '00:15:00'
    mem = 200
    job_bin = '/home/km468/Software/lammps/src/lmp_ufhpc < inp'
    qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                      walltime=walltime,
                                      job_bin=job_bin, mem=mem)
    checkpoint_files = []
    chkpt_file = 'test_lmp.json'
    # setup calibration jobs and run
    cal = CalibrateLammps(parameters, turn_knobs=turn_knobs,
                          qadapter=qadapter, job_cmd = job_cmd,
                          job_dir=job_dir,
                          checkpoint_file=chkpt_file,
                          cal_logger=logger)
    cal.setup()
    cal.run()
    checkpoint_files.append(chkpt_file)
    return checkpoint_files
    

if __name__ == '__main__':
    steps = [test]
    interval = 60
    launch_daemon(steps, interval, ld_logger=logger)
