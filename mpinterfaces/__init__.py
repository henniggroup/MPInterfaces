from __future__ import division, unicode_literals, print_function

__author__ = ", ".join(["Kiran Mathew", "Joshua Gabriel", "Richard G. Hennig"])
__date__ = "Jul 11 2015"
__version__ = "1.1.1"

from pymatgen.matproj.rest import MPRester

def get_struct_from_mp(formula, MAPI_KEY="dwvz2XCFUEI9fJiR", all_structs=False):
    """
    fetches the structure corresponding to the given formula
    from the materialsproject database
    Note: Get the api key from materialsproject website. The one used
    here is nolonger valid.
    Note: The api key can passed to the function or set to the
    environment variable "MAPI_KEY"
    Note: for the given formula there are many structures available, this
    function returns the first one of those structures
    """
    with MPRester(MAPI_KEY) as m:
        data = m.get_data(formula)
        structures = []
        print("\nnumber of structures matching the chemical formula {0} = {1}".format(formula, len(data)) )
        for d in data:
            x = {}
            x['material_id'] = str(d['material_id'])
            structure = m.get_structure_by_material_id(x['material_id'])
            if all_structs:
                structures.append(structure)
            else:
                return structure
        return structures


from .calibrate import Calibrate, CalibrateBulk, CalibrateSlab, CalibrateMolecule
from .instrument import MPINTVaspInputSet, MPINTVaspJob
from .data_processor import MPINTComputedEntry, MPINTVaspDrone, MPINTVasprun
from .interface import Interface, Ligand
from .measurement import Measurement
from .firetasks import MPINTCalibrateTask, MPINTMeasurementTask

