# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

__author__ = ", ".join(
    ["Kiran Mathew", "Joshua J. Gabriel", "Arunima K. Singh", "Michael V. Ashton",\
     "Joshua T. Paul", "Seve G. Moanahan", "Richard G. Hennig"])
__date__ = "March 3 2017"
__version__ = "1.6.0"

import os
import sys
import operator
from pymatgen.matproj.rest import MPRester
import mpinterfaces
from monty.serialization import loadfn

PACKAGE_PATH = mpinterfaces.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')

# set environ variables for MAPI_KEY and VASP_PSP_DIR

try:
    MY_CONFIG = loadfn(PACKAGE_PATH + 'config_mine.yaml')
    try:
        os.environ['VASP_PSP_DIR'] = MY_CONFIG['potentials']
        os.environ['MAPI_KEY'] = MY_CONFIG['mp_api']
    except:
        raise ValueError('config_mine.yaml file not configured .. please'
                         ' set variables potentials and mp_api and retry')

except IOError:
    raise ValueError('No config_mine.yaml file found. Please check')

MAPI_KEY = os.environ.get("MAPI_KEY", "")
MPR = MPRester(MAPI_KEY)

USERNAME= MY_CONFIG['username']
STD_BINARY= MY_CONFIG['normal_binary']
TWOD_BINARY= MY_CONFIG['twod_binary']
VDW_KERNEL= MY_CONFIG['vdw_kernel']
VASP_POTENTIALS= MY_CONFIG['potentials']
QUEUE_SYSTEM= MY_CONFIG['queue_system']

def get_struct_from_mp(formula, MAPI_KEY="", all_structs=False):
    """
    fetches the structure corresponding to the given formula
    from the materialsproject database.

    Note: Get the api key from materialsproject website. The one used
    here is nolonger valid.

    Note: for the given formula there are many structures available,
    this function returns the one with the lowest energy above the hull
    unless all_structs is set to True
    """
    if not MAPI_KEY:
        MAPI_KEY = os.environ.get("MAPI_KEY", "")
        if not MAPI_KEY:
            print('API key not provided')
            print(
                'get API KEY from materialsproject and set it to the MAPI_KEY environment variable. aborting ... ')
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
