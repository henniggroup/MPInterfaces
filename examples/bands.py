# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
reads in KPOINTS(with labels for high symmetry kpoints) and 
vasprun.xml files and plots the band structure along the high 
symmetry kpoints
"""

# To use matplotlib on Hipergator, uncomment the following 2 lines:
# import matplotlib
# matplotlib.use('Agg')

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter

if __name__ == "__main__":
    # readin bandstructure from vasprun.xml and labeled KPOINTS
    run = Vasprun("vasprun.xml", parse_projected_eigen=True)
    bands = run.get_band_structure("KPOINTS", line_mode=True,
                                   efermi=run.efermi)
    bsp = BSPlotter(bands)
    # Blue lines are up spin, red lines are down spin
    bsp.save_plot('band_diagram.eps', ylim=(-5, 5))
    # bsp = BSPlotterProjected(bands)
    # plt = bsp.get_projected_plots_dots( {'Fe':['s', 'p', 'd'],
    #                                     'Sb':['s', 'p', 'd']})
    # get_elt_projected_plots_color()
