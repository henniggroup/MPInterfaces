from __future__ import division, unicode_literals, print_function

import os
import sys
from collections import OrderedDict
import numpy as np

from pymatgen.io.aseio import AseAtomsAdaptor

from mpinterfaces import get_struct_from_mp
from mpinterfaces.lammps import MPINTLammps, CalibrateLammps
from mpinterfaces.utils import get_run_cmmnd


def get_lmp_obj(structure, pair_coeff_file, label):
    """
    return MPINTlammps object
    """
    a, b, c = structure.lattice.abc
    structure.make_supercell([int(12/a), int(12/b), int(12/c)])
    types_of_species = ' '.join( [tos.symbol
                                  for tos in structure.types_of_specie] )
    atomic_mass = [str(i+1)+' '+tos.atomic_mass.__repr__()
                   for i,tos in enumerate(structure.types_of_specie) ]
    parameters = {'atom_style': 'charge',
                  'pair_style' : 'comb3 polar_off',
                  'pair_coeff' : ['* * {0} {1}'.format(pair_coeff_file,
                                                       types_of_species)],
                  'mass': atomic_mass,
                  'minimize':'1.0e-13  1.0e-20  1000  10000',
                  'fix':['fix_nve all nve',
                         '1 all box/relax aniso 0.0 vmax 0.001',
                         '1a all qeq/comb 1 0.0001 file fq.out'] }
    atoms = AseAtomsAdaptor().get_atoms(structure)
    return MPINTLammps(atoms, parameters=parameters, label=label)


# list of structures
structures = get_struct_from_mp('ZnO', all_structs=True)
# path to the coefficient file
pair_coeff_file = os.path.join(os.getcwd(), "ffield.comb3")
#set jobs for the list of MPINTLammps objects
lmp_list = []
for i,structure in enumerate(structures):
    lmp = get_lmp_obj(structure,
                      pair_coeff_file,
                      str(i))
    lmp_list.append(lmp)
turn_knobs = OrderedDict(
    [
        ('LMPS', lmp_list)
    ])
#job directory and run settings
job_dir = 'lammps_job'
qadapter = None
job_cmd = None
nprocs = 8
nnodes = 1
walltime = '00:30:00'
mem = 1000
job_bin = '/home/km468/Software/lammps/src/lmp_ufhpc < inp'
qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                  walltime=walltime,
                                  job_bin=job_bin, mem=mem)
# setup calibration jobs and run
cal = CalibrateLammps(turn_knobs=turn_knobs, qadapter=qadapter,
                      job_cmd = job_cmd, job_dir=job_dir)
cal.setup()
cal.run()
