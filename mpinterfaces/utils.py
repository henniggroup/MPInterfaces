from __future__ import division, unicode_literals, print_function

"""
Utility functions:
Currently creates an ASE slab from a pymatgen input 
structure
"""

import math

from pymatgen.core.surface import Slab
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.aseio import AseAtomsAdaptor

from ase.lattice.surface import surface

def get_ase_slab(pmg_struct, hkl=(1,1,1), min_thick=10, min_vac=10):
    """
    takes in the intial structure as pymatgen Structure object
    uses ase to generate the slab
    returns pymatgen Slab object
    
    Args:
	pmg_struct: pymatgen structure object
	hkl: hkl index of surface of slab to be created
	min_thick: minimum thickness of slab in Angstroms
	min_vac: minimum vacuum spacing 
    """
    ase_atoms = AseAtomsAdaptor().get_atoms(pmg_struct)
    pmg_slab_gen = SlabGenerator(pmg_struct, hkl, min_thick, min_vac)
    h = pmg_slab_gen._proj_height
    nlayers = int(math.ceil(pmg_slab_gen.min_slab_size / h))
    ase_slab = surface(ase_atoms, hkl, nlayers)
    ase_slab.center(vacuum=min_vac/2, axis=2)
    pmg_slab_structure = AseAtomsAdaptor().get_structure(ase_slab)
    return Slab(lattice=pmg_slab_structure.lattice,
                species=pmg_slab_structure.species_and_occu,
                coords=pmg_slab_structure.frac_coords,
                site_properties=pmg_slab_structure.site_properties,
                miller_index=hkl, oriented_unit_cell=pmg_slab_structure,
                shift=0., scale_factor=None, energy=None)
