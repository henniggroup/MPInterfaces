# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
uses openbabel python api to convert a smiles string represntation of
a molecule to its xyz representation

openbabel with python extensions must be made available
"""

try:
    import pybel as pb
except ImportError:
    print("Install openbabel with python bindings and set the path")
    pb = None

smi_input = "C1=CC=CS1"
mol = pb.readstring("smi", smi_input)
mol.make3D()
mol.write(format="xyz", filename='out.xyz', overwrite=True)
