import os
import warnings
from pymatgen.matproj.rest import MPRester

__author__ = "Michael Ashton, Joshua J. Gabriel, Joshua T. Paul, Seve G. Monahan"


try:
    from mpinterfaces import MPINT_CONFIG
    MPR = MPRester(os.environ['MP_API']) if 'MP_API' in os.environ else \
        MPRester(MPINT_CONFIG.get('mp_api', None))
    VASP = MPINT_CONFIG.get('normal_binary', None)
    VASP_2D = MPINT_CONFIG.get('twod_binary', None)
    POTENTIAL_PATH = MPINT_CONFIG.get('potentials')
    USR = MPINT_CONFIG.get('username')
    VDW_KERNEL = MPINT_CONFIG.get('vdw_kernel')
    QUEUE = MPINT_CONFIG.get('queue_system', '').lower()
except ImportError:
    QUEUE = None
    warnings.warn('config_mine.yaml file(required for the twod_materials '
                  'subpackage) not configured. Please set '
                  'variables potentials and mp_api and retry')


if not QUEUE:
    if '/ufrc/' in os.getcwd():
        QUEUE = 'slurm'
    elif '/scratch/' in os.getcwd():
        QUEUE = 'pbs'
    else:
        QUEUE = 'N/A'
