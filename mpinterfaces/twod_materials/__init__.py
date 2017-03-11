import os

from pymatgen.matproj.rest import MPRester

from mpinterfaces import MY_CONFIG

__author__ = "Michael Ashton, Joshua J. Gabriel, " \
             "Joshua T. Paul, Seve G. Monahan"

if 'MP_API' in os.environ:
    MPR = MPRester(os.environ['MP_API'])
else:
    MPR = MPRester(MY_CONFIG['mp_api'])

try:
    VASP = MY_CONFIG['normal_binary']
    VASP_2D = MY_CONFIG['twod_binary']
    POTENTIAL_PATH = MY_CONFIG['potentials']
    USR = MY_CONFIG['username']
    VDW_KERNEL = MY_CONFIG['vdw_kernel']
except:
    raise ValueError('config_mine.yaml file(required for the twod_materials '
                     'subpackage) not configured. Please set '
                     'variables potentials and mp_api and retry')

if 'queue_system' in MY_CONFIG:
    QUEUE = MY_CONFIG['queue_system'].lower()
elif '/ufrc/' in os.getcwd():
    QUEUE = 'slurm'
elif '/scratch/' in os.getcwd():
    QUEUE = 'pbs'
else:
    QUEUE = 'N/A'
