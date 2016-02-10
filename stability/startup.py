import os

import twod_materials.utils as utl

from pymatgen.matproj.rest import MPRester
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints, Incar

from monty.serialization import loadfn

import twod_materials


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')

INCAR_DICT = {
    '@class': 'Incar', '@module': 'pymatgen.io.vasp.inputs', 'AGGAC': 0.0,
    'EDIFF': 1e-06, 'GGA': 'Bo', 'IBRION': 2, 'ISIF': 3, 'ISMEAR': 0,
    'LAECHG': True, 'LCHARG': True, 'LREAL': 'Auto', 'LUSE_VDW': True,
    'NPAR': 4, 'NSW': 50, 'PARAM1': 0.1833333333, 'PARAM2': 0.22,
    'PREC': 'High', 'SIGMA': 0.1, 'LVTOT': True, 'LVHAR': True
    }
KERNEL_PATH = os.path.join(PACKAGE_PATH, 'vdw_kernel.bindat')
MPR = MPRester(
    loadfn(os.path.join(os.path.expanduser('~'), 'config.yaml'))['mp_api']
    )


def relax(directories, submit=True):
    """
    Should be run before pretty much anything else, in order to get the
    right energies of the 2D materials.
    """

    for directory in directories:
        if not utl.is_converged(directory):
            os.chdir(directory)

            # Ensure 20A interlayer vacuum
            utl.add_vacuum(20 - utl.get_spacing(), 0.9)

            # vdw_kernel.bindat file required for VDW calculations.
            os.system('cp {} .'.format(KERNEL_PATH))

            # KPOINTS
            Kpoints.automatic_density(Structure.from_file('POSCAR'),
                                      1000).write_file('KPOINTS')

            # INCAR
            Incar.from_dict(INCAR_DICT).write_file('INCAR')

            # POTCAR
            utl.write_potcar()

            # Submission script
            utl.write_runjob(directory, 1, 8, '600mb', '6:00:00', 'vasp_noz')

            if submit:
                os.system('qsub runjob')
            os.chdir('../')


def relax_competing_species(competing_species, submit=True):
    """
    After obtaining the competing species, relax them with the same
    input parameters as the 2D materials in order to ensure
    compatibility.
    """

    if not os.path.isdir('all_competitors'):
        os.mkdir('all_competitors')
    os.chdir('all_competitors')

    for specie in competing_species:
        if not os.path.isdir(specie[0]):
            os.mkdir(specie[0])
        if not utl.is_converged(specie[0]):
            os.chdir(specie[0])
            os.system('cp {} .'.format(KERNEL_PATH))
            structure = MPR.get_structure_by_material_id(specie[1])
            structure.to('POSCAR', 'POSCAR')
            Kpoints.automatic_density(structure, 1000).write_file('KPOINTS')
            Incar.from_dict(INCAR_DICT).write_file('INCAR')
            utl.write_potcar()
            utl.write_runjob(specie[0], 1, 8, '600mb', '6:00:00', 'vasp')
            if submit:
                os.system('qsub runjob')
            os.chdir('../')
    os.chdir('../')
