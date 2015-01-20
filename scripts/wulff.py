"""
Construct nanoparticle given surface miller indices and the corresponding surface energies
uses wulff construction
requires the ase package
"""

from ase.cluster import wulff_construction
from pymatgen.io.aseio import AseAtomsAdaptor

symbol = 'Pt'
surfaces = [ (1,0,0), (1,1,1) ]
surface_energies = [1, 1]
size = 200 #number of atoms
structure = "fcc"
latticeconstant = 5.0
atoms = wulff_construction(symbol, surfaces, surface_energies, size, structure,
                           rounding='closest', latticeconstant=latticeconstant,
                           debug=False, maxiter=100)
#convert to pymatgen structure
pgen_structure = AseAtomsAdaptor().get_structure(atoms)
pgen_structure.to(fmt='poscar', filename='POSCAR_pt_nano.vasp')
