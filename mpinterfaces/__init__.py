from __future__ import division, unicode_literals, print_function

__author__ = ", ".join(["Kiran Mathew", "Joshua Gabriel", "Arunima Singh", "Richard G. Hennig"])
__date__ = "Aug 10 2015"
__version__ = "1.1.2"

import os
import sys
import operator
from pymatgen.matproj.rest import MPRester

def get_struct_from_mp(formula, MAPI_KEY="", all_structs=False):
    """
    fetches the structure corresponding to the given formula
    from the materialsproject database.
    
    Note: Get the api key from materialsproject website. The one used
    here is nolonger valid.
    
    Note: for the given formula there are many structures available, this
    function returns the one with the lowest energy above the hull
    unless all_structs is set to True
    """
    if not MAPI_KEY:
        MAPI_KEY = os.environ.get("MAPI_KEY", "")
        if not MAPI_KEY:
            print('API key not provided')
            print('get API KEY from materialsproject and set it to the MAPI_KEY environment variable. aborting ... ')
            sys.exit()                            
    with MPRester(MAPI_KEY) as m:
        data = m.get_data(formula)
        structures = []
        x = {}
        print("\nnumber of structures matching the chemical formula {0} = {1}".format(formula, len(data)) )
        print("The one with the the lowest energy above the hull is returned, unless all_structs is set to True")
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
            print("The id of the material corresponding to the lowest energy above the hull = {0}".format(mineah_key))
            if mineah_key:
                return m.get_structure_by_material_id(mineah_key)
            else:
                return None


from .calibrate import Calibrate, CalibrateBulk, CalibrateSlab, CalibrateMolecule
from .instrument import MPINTVaspInputSet, MPINTVaspJob
from .data_processor import MPINTComputedEntry, MPINTVaspDrone, MPINTVasprun
from .interface import Interface, Ligand
from .measurement import Measurement
from .firetasks import MPINTCalibrateTask, MPINTMeasurementTask

