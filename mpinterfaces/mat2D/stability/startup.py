from __future__ import print_function, division, unicode_literals

import os

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints, Incar

from mpinterfaces import VASP_TWOD_BIN, VDW_KERNEL, QUEUE_SYSTEM
import mpinterfaces.utils as utl

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

INCAR_DICT = {
    '@class': 'Incar', '@module': 'pymatgen.io.vasp.inputs', 'AGGAC': 0.0,
    'EDIFF': 1e-04, 'GGA': 'Bo', 'IBRION': 2, 'ISIF': 3, 'ISMEAR': 1,
    'LAECHG': True, 'LCHARG': True, 'LREAL': 'Auto', 'LUSE_VDW': True,
    'NPAR': 4, 'NSW': 50, 'PARAM1': 0.1833333333, 'PARAM2': 0.22,
    'PREC': 'Accurate', 'ENCUT': 500, 'SIGMA': 0.1, 'LVTOT': True,
    'LVHAR': True, 'ALGO': 'Fast', 'ISPIN': 2
    }


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
        if VDW_KERNEL:
            os.system('cp {} .'.format(VDW_KERNEL))
        # KPOINTS
        Kpoints.automatic_density(Structure.from_file('POSCAR'),
                                  1000).write_file('KPOINTS')

        # INCAR
        INCAR_DICT.update(
            {'MAGMOM': utl.get_magmom_string(Structure.from_file('POSCAR'))}
        )
        Incar.from_dict(INCAR_DICT).write_file('INCAR')
        # POTCAR
        utl.write_potcar()

        # Special tasks only performed for 2D materials.
        if dim == 2:
            # Ensure 20A interlayer vacuum
            utl.ensure_vacuum(Structure.from_file('POSCAR'), 20)
            # Remove all z k-points.
            kpts_lines = open('KPOINTS').readlines()
            with open('KPOINTS', 'w') as kpts:
                for line in kpts_lines[:3]:
                    kpts.write(line)
                kpts.write(kpts_lines[3].split()[0] + ' '
                           + kpts_lines[3].split()[1] + ' 1')

        # Submission script
        if QUEUE_SYSTEM == 'pbs':
            utl.write_pbs_runjob(directory, 1, 16, '800mb', '6:00:00',
                VASP_TWOD_BIN)
            submission_command = 'qsub runjob'

        elif QUEUE_SYSTEM == 'slurm':
            utl.write_slurm_runjob(directory, 16, '800mb', '6:00:00',
                VASP_TWOD_BIN)
            submission_command = 'sbatch runjob'

        if submit:
            os.system(submission_command)
