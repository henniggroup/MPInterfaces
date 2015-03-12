from __future__ import division, unicode_literals, print_function

"""
uses openbabel python api to convert a smiles string represntation of
a molecule to its xyz representation

openbabel with python extensions must be made available
"""

import openbabel
import pybel
from pybel import readstring

smi_input = "C1=CC=CS1"
mol = readstring("smi", smi_input) 
mol.make3D()
mol.write(format="xyz", filename='out.xyz', overwrite=True)
