import os

from pymatgen.io.vasp.inputs import Incar

from mpinterfaces.twod_materials.stability import get_magmom_string, INCAR_DICT

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__version__ = "1.6"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


def run_major_axis_anisotropy_calculations(submit=True):
    """
    Perform static calculations with the magnetic axis along
    100, 010, and 001.

    Kwargs:
        submit (bool): Whether or not to submit the job.
    """

    if not os.path.isdir('MAE'):
        os.mkdir('MAE')
    os.chdir('MAE')
    for d in ['100', '010', '001']:
        if not os.path.isdir(d):
            os.path.mkdir(d)
        os.chdir(d)
        os.system('cp ../CONTCAR POSCAR')
        os.system('cp ../POTCAR .')
        axis = [float(char) for char in d]
        # Small positive number, see vasp manual
        if d in ['001', '010']:
            axis[0] = 0.00000001
        else:
            axis[1] = 0.00000001

        saxis = ' '.join(axis)
        incar_dict = INCAR_DICT
        incar_dict.update({'EDIFF': 1e-8, 'GGA_COMPAT': False, 'ISMEAR': -5,
                           'LORBIT': 11, 'LSORBIT': True, 'LWAVE': False,
                           'LCHARG': False, 'LAECHG': False,
                           'MAGMOM': get_magmom_string(), 'SAXIS': saxis})
        Incar.from_dict(incar_dict).write_file('INCAR')

        if submit:
            os.system(submission_command)


def run_xy_anisotropy_calculations(resolution=10, submit=True):
    """
    Perform static calculations with the magnetic axis along
    several directions between 100 and 010.

    Kwargs:
        resolution (int): step size between axes. The total
            number of calculations will be 90 / `resolution`.
        submit (bool): Whether or not to submit the job.
    """


def run_xz_anisotropy_calculations(resolution=10, submit=True):
    """
    Perform static calculations with the magnetic axis along
    several directions between 100 and 001.

    Kwargs:
        resolution (int): step size between axes. The total
            number of calculations will be 90 / `resolution`.
        submit (bool): Whether or not to submit the job.
    """
