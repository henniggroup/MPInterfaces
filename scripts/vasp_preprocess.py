"""
This script demonstrates the usage of the module mpinterfaces/instrument.py

Creates a sample workflow for ENCUT convergence study
"""

import os
import sys
import shutil as shu
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from pymatgen.core.structure import Structure
import subprocess as sp
from mpinterfaces.instrument import CalibrateMolecule

#dictionary of defaults for molecule ENCUT convergence
default_incar = {'SYMPREC': 1e-06, 
                 'SYSTEM': 'molecule', 
                 'ENCUT': 500, 
                 'ISIF': 2, 
                 'IBRION': 2, 
                 'ALGO': 'Normal', 
                 'ISMEAR': 1, 
                 'ISPIN': 1, 
                 'EDIFF': 1e-06, 
                 'EDIFFG': -0.001, 
                 'NPAR': 8, 
                 'SIGMA': 0.1, 
                 'PREC': 'Accurate'}

incar = Incar.from_dict(default_incar)
poscar = Poscar.from_file('POSCAR', check_for_POTCAR=False)#Structure.from_file("POSCAR", sort=True)
symbols = poscar.site_symbols #symbol_set
potcar = Potcar(symbols)
kpoints = Kpoints()

calmol = CalibrateMolecule(incar, poscar, potcar, kpoints)
calmol.encut_cnvg(range(400,800,100))
calmol.run()
