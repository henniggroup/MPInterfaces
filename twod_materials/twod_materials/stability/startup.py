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

INCAR_DICT = {
    '@class': 'Incar', '@module': 'pymatgen.io.vasp.inputs', 'AGGAC': 0.0,
    'EDIFF': 1e-06, 'GGA': 'Bo', 'IBRION': 2, 'ISIF': 3, 'ISMEAR': 1,
    'LAECHG': True, 'LCHARG': True, 'LREAL': 'Auto', 'LUSE_VDW': True,
    'NPAR': 4, 'NSW': 50, 'PARAM1': 0.1833333333, 'PARAM2': 0.22,
    'PREC': 'Accurate', 'ENCUT': 500, 'SIGMA': 0.1, 'LVTOT': True,
    'LVHAR': True, 'ALGO': 'Fast', 'ISPIN': 2
    }
KERNEL_PATH = os.path.join(PACKAGE_PATH, 'vdw_kernel.bindat')

if '/ufrc/' in os.getcwd():
    HIPERGATOR = 2
elif '/scratch/' in os.getcwd():
    HIPERGATOR = 1

try:
    MPR = MPRester(
        loadfn(os.path.join(os.path.expanduser('~'), 'config.yaml'))['mp_api']
        )
    VASP = loadfn(os.path.join(os.path.expanduser('~'),
                               'config.yaml'))['normal_binary']
    VASP_2D = loadfn(os.path.join(os.path.expanduser('~'),
                                  'config.yaml'))['twod_binary']
except IOError:
    try:
        MPR = MPRester(
            os.environ['MP_API']
            )
    except KeyError:
        raise ValueError('No config.yaml file found. Please check'
                         ' that your config.yaml is in your home directory'
                         ' and contains the field'
                         ' mp_api: your_api_key')


def get_magmom_string():
    """
    Based on a POSCAR, returns the string required for the MAGMOM
    setting in the INCAR. Initializes transition metals with 6.0
    bohr magneton and all others with 0.5.
    """

    magmoms = []
    poscar_lines = open('POSCAR').readlines()
    elements = poscar_lines[5].split()
    amounts = poscar_lines[6].split()
    for i in range(len(elements)):
        if Element(elements[i]).is_transition_metal:
            magmoms.append('{}*6.0'.format(amounts[i]))
        else:
            magmoms.append('{}*0.5'.format(amounts[i]))
    return ' '.join(magmoms)


def relax(submit=True, force_overwrite=False, dim='2D'):
    """
    Should be run before pretty much anything else, in order
    to get the right energy of the 2D material.

    `dim`: Set to "3D" if you want to include 3D k-points. The
    default behavior is for 2D materials, and removes k-points with
    z-components.
    """

    if force_overwrite or not utl.is_converged(os.getcwd()):
        directory = os.getcwd().split('/')[-1]
        # Ensure 20A interlayer vacuum
        utl.add_vacuum(20 - utl.get_spacing(), 0.9)
        # vdw_kernel.bindat file required for VDW calculations.
        os.system('cp {} .'.format(KERNEL_PATH))
        # KPOINTS
        Kpoints.automatic_density(Structure.from_file('POSCAR'),
                                  1000).write_file('KPOINTS')
        if dim == '2D':
            kpts_lines = open('KPOINTS').readlines()
            with open('KPOINTS', 'w') as kpts:
                for line in kpts_lines[:3]:
                    kpts.write(line)
                kpts.write(kpts_lines[3].split()[0] + ' '
                           + kpts_lines[3].split()[1] + ' 1')
        # INCAR
        INCAR_DICT.update({'MAGMOM': get_magmom_string()})
        Incar.from_dict(INCAR_DICT).write_file('INCAR')
        # POTCAR
        utl.write_potcar()
        # Submission script
        if HIPERGATOR == 1:
            utl.write_pbs_runjob(directory, 1, 16, '800mb', '6:00:00',
                                 VASP_2D)
            submission_command = 'qsub runjob'

        elif HIPERGATOR == 2:
            utl.write_slurm_runjob(directory, 16, '800mb', '6:00:00',
                                   VASP_2D)
            submission_command = 'sbatch runjob'

        if submit:
            os.system(submission_command)


def relax_competing_phases(competing_phases, submit=True,
                            force_overwrite=False):
    """
    After obtaining the competing phases, relax them with the same
    input parameters as the 2D materials in order to ensure
    compatibility.
    """

    if not os.path.isdir('all_competitors'):
        os.mkdir('all_competitors')
    os.chdir('all_competitors')

    for phase in competing_phases:
        if not os.path.isdir(phase[0]):
            os.mkdir(phase[0])
        directory = os.path.join(os.getcwd(), phase[0])
        if force_overwrite or not utl.is_converged(directory):
            os.chdir(phase[0])
            os.system('cp {} .'.format(KERNEL_PATH))
            structure = MPR.get_structure_by_material_id(phase[1])
            structure.to('POSCAR', 'POSCAR')
            Kpoints.automatic_density(structure, 1000).write_file('KPOINTS')
            INCAR_DICT.update({'MAGMOM': get_magmom_string()})
            Incar.from_dict(INCAR_DICT).write_file('INCAR')
            utl.write_potcar()
            if HIPERGATOR == 1:
                utl.write_pbs_runjob('{}_3d'.format(phase[0]), 1, 8, '600mb',
                                     '6:00:00', VASP)
                submission_command = 'qsub runjob'

            elif HIPERGATOR == 2:
                utl.write_slurm_runjob('{}_3d'.format(phase[0]), 8, '600mb',
                                       '6:00:00', VASP)
                submission_command = 'sbatch runjob'

            if submit:
                os.system(submission_command)

            os.chdir('../')
    os.chdir('../')


def relax_3d(submit=True, force_overwrite=False):
    """
    Standard relaxation for a single directory of a bulk material.
    """

    if force_overwrite or not utl.is_converged(os.getcwd()):
        directory = os.getcwd().split('/')[-1]

        # vdw_kernel.bindat file required for VDW calculations.
        os.system('cp {} .'.format(KERNEL_PATH))
        # KPOINTS
        Kpoints.automatic_density(Structure.from_file('POSCAR'),
                                  1000).write_file('KPOINTS')
        # INCAR
        INCAR_DICT.update({'MAGMOM': get_magmom_string()})
        Incar.from_dict(INCAR_DICT).write_file('INCAR')
        # POTCAR
        utl.write_potcar()
        # Submission script
        if HIPERGATOR == 1:
            utl.write_pbs_runjob('{}_3d'.format(directory), 1, 8, '600mb',
                                 '6:00:00', VASP)
            submission_command = 'qsub runjob'

        elif HIPERGATOR == 2:
            utl.write_slurm_runjob('{}_3d'.format(directory), 8, '600mb',
                                   '6:00:00', VASP)
            submission_command = 'sbatch runjob'

        if submit:
            os.system(submission_command)
