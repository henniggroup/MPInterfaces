# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

""" use the interface module to create bare slab """

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar

from mpinterfaces.interface import Interface

if __name__ == '__main__':
    # input structure file: bulk
    fin = 'POSCAR.bulk'
    # output file name
    fout = 'POSCAR'
    hkl = [0, 0, 1]  # hkl wrt the input structure
    min_thick = 15  # angstroms
    min_vac = 30  # angstroms
    # use ase to create an orthogonal slab
    bulk = Structure.from_file(fin)
    iface = Interface(bulk, hkl=hkl,
                      min_thick=min_thick, min_vac=min_vac,
                      primitive=False, from_ase=True)
    iface.create_interface()
    iface.sort()
    # set selective dynamics flags
    # 1 --> T and 0 --> F
    sd_flags = [[1, 1, 1] for i in iface.sites]
    iface_poscar = Poscar(iface, selective_dynamics=sd_flags)
    # write to file
    iface_poscar.write_file(fout)
