from pymatgen.matproj.rest import MPRester

def get_struct_from_mp(formula):
    """
    fetches the structure corresponding to the given formula
    from the materialsproject database
    Note: get the api key from materialsproject website
    provide the api key here os set the environment variable "MAPI_KEY"
    Note: for the given formula there are many structures available, this
    function returns the first one of those structures
    """
    with MPRester("dwvz2XCFUEI9fJiR") as m:
        data = m.get_data(formula)
        print "\nnumber of structures matching the chemical formula "+formula+" = ", len(data)
        for d in data:
            x = {}
            x['material_id'] = str(d['material_id'])
            structure = m.get_structure_by_material_id(x['material_id'])
            return structure


from .calibrate import Calibrate, CalibrateBulk, CalibrateSlab, CalibrateMolecule
from .instrument import MPINTVaspInputSet, MPINTVaspJob
from .data_processor import MPINTComputedEntry, MPINTVaspDrone, MPINTVasprun
from .interface import Interface, Ligand
from .measurement import Measurement
from .firetasks import MPINTCalibrateTask, MPINTMeasurementTask

