# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Sample LAMMPS workflow that runs MD calculations for a bunch of structures and 
subsequently generate their phasediagram.
"""

from math import ceil

import matplotlib

matplotlib.use('Agg')

from pymatgen.core.composition import Composition
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter

from mpinterfaces import get_struct_from_mp
from mpinterfaces.lammps import CalibrateLammps
from mpinterfaces.utils import *

# all the info/warnings/outputs redirected to the log file: 
# lammps_Al_O.log
logger = get_logger('lammps_Al_O')
# list of structures from materialsproject
structures = get_struct_from_mp('Al-O', all_structs=True)
# scale the structures
scell_size = 12
for s in structures:
    a, b, c = s.lattice.abc
    s.make_supercell([ceil(scell_size / a),
                      ceil(scell_size / b),
                      ceil(scell_size / c)])
# lammps input paramaters    
parameters = {'atom_style': 'charge',
              'charges': {'Al': 0, 'O': 0},
              'minimize': '1.0e-13  1.0e-20  1000  10000',
              'fix': ['fix_nve all nve',
                      '1 all box/relax aniso 0.0 vmax 0.001',
                      '1a all qeq/comb 1 0.0001 file fq.out']}
# list of pair styles
pair_styles = ['comb3 polar_off']
# list of pair coefficient files
# this file must be in the folder where this script is run
pair_coeff_files = [os.path.join(os.getcwd(), "ffield.comb3")]


def step1(**kwargs):
    """
    setup and run all lammps jobs
    """
    turn_knobs = OrderedDict(
        [
            ('STRUCTURES', structures),
            ('PAIR_STYLE', pair_styles),
            ('PAIR_COEFF', pair_coeff_files)
        ])
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
    chkpt_file = 'step1.json'
    # setup calibration jobs and run
    cal = CalibrateLammps(parameters, turn_knobs=turn_knobs,
                          qadapter=qadapter, job_cmd=job_cmd,
                          job_dir=job_dir, is_matrix=True,
                          checkpoint_file=chkpt_file,
                          cal_logger=logger)
    cal.setup()
    cal.run()
    checkpoint_files.append(chkpt_file)
    return checkpoint_files


def step2(**kwargs):
    """
    post process:
       get energies from the jobs in the previous step and 
       generate the phase diagram
    """
    chkfile = kwargs['checkpoint_files'][0]
    all_jobs = jobs_from_file(chkfile)
    entries = []
    # add endpoint data
    Al = Composition("Al1O0")
    energy_al = -3.36
    O = Composition("Al0O1")
    energy_o = -2.58
    entries.append(PDEntry(Al, energy_al))
    entries.append(PDEntry(O, energy_o))
    # get data and create entries
    for job in all_jobs:
        comp = job.vis.mplmp.structure.composition
        energy = job.final_energy
        entries.append(PDEntry(comp, energy))
    pd = PhaseDiagram(entries)
    plotter = PDPlotter(pd, show_unstable=True)
    plotter.write_image('Al_O_phasediagram.jpg')
    return None


if __name__ == '__main__':
    steps = [step1, step2]
    interval = 60
    launch_daemon(steps, interval, ld_logger=logger)
