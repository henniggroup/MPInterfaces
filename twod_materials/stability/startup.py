from __future__ import print_function, division, unicode_literals

import os

import twod_materials.utils as utl

from pymatgen.matproj.rest import MPRester
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Kpoints, Incar

from monty.serialization import loadfn

import twod_materials


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')
PACKAGE_PATH = '/'.join(PACKAGE_PATH.split('/')[:-2])

INCAR_DICT = {
    '@class': 'Incar', '@module': 'pymatgen.io.vasp.inputs', 'AGGAC': 0.0,
    'EDIFF': 1e-04, 'GGA': 'Bo', 'IBRION': 2, 'ISIF': 3, 'ISMEAR': 1,
    'LAECHG': True, 'LCHARG': True, 'LREAL': 'Auto', 'LUSE_VDW': True,
    'NPAR': 4, 'NSW': 50, 'PARAM1': 0.1833333333, 'PARAM2': 0.22,
    'PREC': 'Accurate', 'ENCUT': 500, 'SIGMA': 0.1, 'LVTOT': True,
    'LVHAR': True, 'ALGO': 'Fast', 'ISPIN': 2
    }

try:
    config_vars = loadfn(os.path.join(os.path.expanduser('~'), 'config.yaml'))
except:
    print('WARNING: No config.yaml file was found. Please configure the '\
    'config.yaml and put it in your home directory.')
    # Still set them for testing purposes.
    config_vars = loadfn(os.path.join(PACKAGE_PATH, 'config.yaml'))
if 'MP_API' in os.environ:
    MPR = MPRester(os.environ['MP_API'])
else:
    MPR = MPRester(config_vars['mp_api'])
VASP = config_vars['normal_binary']
VASP_2D = config_vars['twod_binary']
VDW_KERNEL = config_vars['vdw_kernel']
if 'queue_system' in config_vars:
    QUEUE = config_vars['queue_system'].lower()
elif '/ufrc/' in os.getcwd():
    QUEUE = 'slurm'
elif '/scratch/' in os.getcwd():
    QUEUE = 'pbs'


def relax(dim=2, submit=True, force_overwrite=False):
    """
    Writes input files and (optionally) submits a self-consistent
    relaxation. Should be run before pretty much anything else, in
    order to get the right energy and structure of the material.

    Args:
        dim (int): 2 for relaxing a 2D material, 3 for a 3D material.
        submit (bool): Whether or not to submit the job.
        force_overwrite (bool): Whether or not to overwrite files
            if an already converged vasprun.xml exists in the
            directory.
    """

    if force_overwrite or not utl.is_converged(os.getcwd()):
        directory = os.getcwd().split('/')[-1]

        # vdw_kernel.bindat file required for VDW calculations.
        if VDW_KERNEL != '/path/to/vdw_kernel.bindat':
            os.system('cp {} .'.format(VDW_KERNEL))
        # KPOINTS
        Kpoints.automatic_density(Structure.from_file('POSCAR'),
                                  1000).write_file('KPOINTS')

        # INCAR
        INCAR_DICT.update({'MAGMOM': utl.get_magmom_string()})
        Incar.from_dict(INCAR_DICT).write_file('INCAR')
        # POTCAR
        utl.write_potcar()

        # Special tasks only performed for 2D materials.
        if dim == 2:
            # Ensure 20A interlayer vacuum
            utl.add_vacuum(20 - utl.get_spacing(), 0.9)
            # Remove all z k-points.
            kpts_lines = open('KPOINTS').readlines()
            with open('KPOINTS', 'w') as kpts:
                for line in kpts_lines[:3]:
                    kpts.write(line)
                kpts.write(kpts_lines[3].split()[0] + ' '
                           + kpts_lines[3].split()[1] + ' 1')

        # Submission script
        if QUEUE == 'pbs':
            utl.write_pbs_runjob(directory, 1, 16, '800mb', '6:00:00',
                                 VASP_2D)
            submission_command = 'qsub runjob'

        elif QUEUE == 'slurm':
            utl.write_slurm_runjob(directory, 16, '800mb', '6:00:00',
                                   VASP_2D)
            submission_command = 'sbatch runjob'

        if submit:
            os.system(submission_command)
