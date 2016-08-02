# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Workflow to study the convergence of VASP parameters: ENCUT and KPOINTS
for a range of systems.
"""

from six.moves import range

from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.inputs import Kpoints

from mpinterfaces import get_struct_from_mp
from mpinterfaces.calibrate import Calibrate
from mpinterfaces.utils import *

# all the info/warnings/outputs redirected to the log file: convg.log
logger = get_logger('convg')

incar_dict = dict(
    PREC='Accurate',
    ENCUT=400,
    ISMEAR=1,
    EDIFF='1E-6',
    NSW=0,
    NPAR=4,
    LCHARG='.FALSE.',
    LWAVE='.FALSE.')
# INCAR
incar = Incar.from_dict(incar_dict)
# KPOINTS
kpoints = Kpoints.monkhorst_automatic(kpts=(12, 12, 12))
# QUE
nprocs = 8
nnodes = 1
mem = '1000'
walltime = '1:00:00'
job_bin = '/home/km468/Software/VASP/vasp.5.3.5/vasp'
qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                  walltime=walltime,
                                  job_bin=job_bin, mem=mem)
# STRUCTURES
structures = ['Pt', 'Ag', 'Cu']
poscar_list = []


def get_structures():
    for species in structures:
        struct_species = get_struct_from_mp(species)
        poscar_list.append(Poscar(struct_species))
    return None


# step 1
def convergence(**kwargs):
    """
    Set up and run all required jobs for the convergence study.
    Recursive(is_matrix=True) ENCUT and KPOINTS jobs for the given list of
    POSCARS.
    Note: Since this is a convergence study ENCUT jobs and KPOINTS jobs are
    independently i.e recursion is used for the list of POSCARS+ENCUT and
    POSCARS+KPOINTS separately
    """
    job_dir = "convergence"
    kpoints_list = [[x, x, x] for x in range(8, 15)]
    encut_list = range(400, 800, 100)
    params = {'ENCUT': encut_list, 'KPOINTS': kpoints_list}
    checkpoint_files = []
    for k, v in params.items():
        turn_knobs = OrderedDict([('POSCAR', poscar_list)])
        turn_knobs[k] = v
        chkpt_file = k + '.json'
        run_cal(turn_knobs, qadapter, job_cmd, job_dir,
                chkpt_file, incar=incar, kpoints=kpoints, is_matrix=True)
        checkpoint_files.append(chkpt_file)
    return checkpoint_files


# step 2
def post_process(**kwargs):
    """
    process the finished jobs from the previous step and obtain the optimum
    value for the parameters.
    """
    # parameters that need to be optimized
    params = ['ENCUT', 'KPOINTS']
    for chkfile in kwargs['checkpoint_files']:
        logger.info('processing {}'.format(chkfile))
        conv_data = get_convergence_data(chkfile, params=params)
        for k in conv_data.keys():
            for param in params:
                optimum = get_opt_params(conv_data, k, param=param,
                                         ev_per_atom=0.001)
                logger.info('For {0}, the  optimum {1} : {2}'.format(k,
                                                                     param,
                                                                     optimum))
    return None


def run_cal(turn_knobs, qadapter, job_cmd, job_dir, checkpoint_file,
            incar=None, poscar=None, potcar=None, kpoints=None,
            Grid_type='G', functional='PBE', is_matrix=True):
    """
    setup and run calibrate job
    """
    cal = Calibrate(incar, poscar, potcar, kpoints,
                    is_matrix=is_matrix,
                    turn_knobs=turn_knobs, qadapter=qadapter,
                    job_cmd=job_cmd, job_dir=job_dir,
                    Grid_type=Grid_type, functional=functional,
                    checkpoint_file=checkpoint_file, cal_logger=logger)
    cal.setup()
    cal.run()


if __name__ == '__main__':
    get_structures()
    steps = [convergence, post_process]
    interval = 30
    launch_daemon(steps, interval, ld_logger=logger)
