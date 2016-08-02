# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
reads in vasprun.xml file and plots the density of states 
"""

# To use matplotlib on Hipergator, uncomment the following 2 lines:
# import matplotlib
# matplotlib.use('Agg')

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter

if __name__ == "__main__":
    # readin the density of states from vasprun.xml file
    run = Vasprun("vasprun.xml", parse_projected_eigen=True)
    complete_dos = run.complete_dos
    print('cbm and vbm ', complete_dos.get_cbm_vbm())
    print('gap = ', complete_dos.get_gap())
    # get orbital projected DOS.    
    spd_dos = complete_dos.get_spd_dos()
    plotter = DosPlotter()
    plotter.add_dos_dict(spd_dos)
    plotter.save_plot('dos.eps')
