# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import os
import sys
import operator
import warnings

__author__ = "Kiran Mathew, Joshua J. Gabriel, Michael Ashton, " \
             "Arunima K. Singh, Joshua T. Paul, Seve G. Monahan, " \
             "Richard G. Hennig"
__date__ = "March 3 2017"
__version__ = "1.7.0"

PACKAGE_PATH = os.path.dirname(os.path.abspath(__file__))

from .config_loader import CONFIG
from pymatgen.ext.matproj import MPRester

# First, correctly fetch PMG_MAPI_KEY using .get() to handle cases where it might not be set
PMG_MAPI_KEY = CONFIG.get('PMG_MAPI_KEY')

MPR = MPRester(PMG_MAPI_KEY)

# Other configurations
USERNAME = CONFIG.get('username')
VASP_STD_BIN = CONFIG.get('normal_binary')
VASP_TWOD_BIN = CONFIG.get('twod_binary')
VDW_KERNEL = CONFIG.get('vdw_kernel')
VASP_PSP = CONFIG.get('potentials')
QUEUE_SYSTEM = CONFIG.get('queue_system')
QUEUE_TEMPLATE = CONFIG.get('queue_template')

if not QUEUE_SYSTEM:
    QUEUE_SYSTEM = 'slurm'
    
def get_struct_from_mp(formula, PMG_MAPI_KEY="", all_structs=False):
    """
    fetches the structure corresponding to the given formula
    from the materialsproject database.

    Note: Get the api key from materialsproject website. The one used
    here is nolonger valid.

    Note: for the given formula there are many structures available,
    this function returns the one with the lowest energy above the hull
    unless all_structs is set to True
    """
    if not PMG_MAPI_KEY:
        PMG_MAPI_KEY = os.environ.get("PMG_MAPI_KEY", "")
        if not PMG_MAPI_KEY:
            print('API key not provided')
            print(
                'get API KEY from materialsproject and set it to the PMG_MAPI_KEY environment variable. aborting ... ')
            sys.exit()
    with MPR as m:
        data = m.get_data(formula)
        structures = []
        x = {}
        print(
            "\nnumber of structures matching the chemical formula {0} = {1}".format(
                formula, len(data)))
        print(
            "The one with the the lowest energy above the hull is returned, unless all_structs is set to True")
        for d in data:
            mpid = str(d['material_id'])
            x[mpid] = d['e_above_hull']
            if all_structs:
                structure = m.get_structure_by_material_id(mpid)
                structures.append(structure)
        if all_structs:
            return structures
        else:
            mineah_key = sorted(x.items(), key=operator.itemgetter(1))[0][0]
            print(
                "The id of the material corresponding to the lowest energy above the hull = {0}".format(
                    mineah_key))
            if mineah_key:
                return m.get_structure_by_material_id(mineah_key)
            else:
                return None
