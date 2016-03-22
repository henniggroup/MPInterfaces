# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
reads in KPOINTS(with labels) and vasprun.xml files and
plots band diagram and density of states 
from http://gvallver.perso.univ-pau.fr/?p=587
"""

from six.moves import range
from six.moves import zip

import numpy as np

try:
    # To use matplotlib on Hipergator, uncomment the following 2 lines:
    # import matplotlib
    # matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.gridspec import GridSpec
except ImportError:
    print("Install matplotlib")
    plt = None

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin


def rgbline(ax, k, e, red, green, blue, alpha=1.):
    # creation of segments based on
    # http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    # amount of r,g,b weighted by the s, p, d contribution to the density
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

    nseg = len(k) - 1
    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=2)
    ax.add_collection(lc)


if __name__ == "__main__":
    # --------------------------------------------------------------  
    # read data
    # --------------------------------------------------------------      
    # readin bandstructure and density of states from vasprun.xml file
    run = Vasprun("vasprun.xml", parse_projected_eigen=True)
    bands = run.get_band_structure("KPOINTS", line_mode=True,
                                   efermi=run.efermi)
    complete_dos = run.complete_dos
    print('cbm and vbm ', complete_dos.get_cbm_vbm())
    print('gap = ', complete_dos.get_gap())
    # get orbital projected DOS.    
    spd_dos = complete_dos.get_spd_dos()
    # kpoints labels, must conform with the label in the KPOINTS file
    labels = [r"$G$", r"$X$", r"$M$", r"$G$"]

    # --------------------------------------------------------------  
    # compute a dictionary of projections on elements and specific orbitals
    # A dictionary of Elements and Orbitals for which we want
    # to have projections on. It is given as: {Element:[orbitals]},
    # e.g., {'Cu':['d','s']}
    # --------------------------------------------------------------          
    pbands = bands.get_projections_on_elts_and_orbitals(
        {"Fe": ["s", "p", "d"]})
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    # loop over the energies and kpoints
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            sc = pbands[Spin.up][b][k]["Fe"]["s"] ** 2
            pc = pbands[Spin.up][b][k]["Fe"]["p"] ** 2
            dc = pbands[Spin.up][b][k]["Fe"]["d"] ** 2
            tot = sc + pc + dc
            if tot != 0.0:
                contrib[b, k, 0] = sc / tot  # normalized s contribution
                contrib[b, k, 1] = pc / tot  # normalized p contribution
                contrib[b, k, 2] = dc / tot  # normalized d contribution

    # --------------------------------------------------------------  
    # set up matplotlib plot
    # --------------------------------------------------------------  
    # general options for plot
    font = {'family': 'serif', 'size': 24}
    plt.rc('font', **font)

    # set up 2 graph with aspec ration 2/1
    # plot 1: bands diagram
    # plot 2: Density of States
    gs = GridSpec(1, 2, width_ratios=[2, 1])
    fig = plt.figure(figsize=(11.69, 8.27))
    fig.suptitle("HSE Band Structure of FeS-litharge")
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])  # , sharey=ax1)

    # --------------------------------------------------------------
    # energies on the y axis    
    # set ylimit for the plot
    # max and min energy values computed from the read in data
    # -------------------------------------------------------------- 
    emin = 1e100
    emax = -1e100
    for spin in bands.bands.keys():
        for b in range(bands.nb_bands):
            emin = min(emin, min(bands.bands[spin][b]))
            emax = max(emax, max(bands.bands[spin][b]))
    emin -= bands.efermi + 1
    emax -= bands.efermi - 1
    emin1 = -5
    emax1 = 5
    ax1.set_ylim(emin1, emax1)
    ax2.set_ylim(emin, emax)

    # -------------------------------------------------------------- 
    # Band Diagram
    # --------------------------------------------------------------
    # plot bands using rgb mapping
    for b in range(bands.nb_bands):
        rgbline(ax1, list(range(len(bands.kpoints))),
                [e - bands.efermi for e in bands.bands[Spin.up][b]],
                contrib[b, :, 0], contrib[b, :, 1], contrib[b, :, 2])
        rgbline(ax1, list(range(len(bands.kpoints))),
                [e - bands.efermi for e in bands.bands[Spin.down][b]],
                contrib[b, :, 0], contrib[b, :, 1], contrib[b, :, 2],
                alpha=0.2)

    # ax1.set_title("Band diagram")       
    ax1.set_xlim(0, len(bands.kpoints))
    ax1.set_xlabel("k-points")
    ax1.set_ylabel(r"$E - E_f$ / eV")
    ax1.grid()
    # fermi level at 0, draw a line through it
    ax1.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="k", lw=2)
    # kpoint labels
    # draw vertical lines along each kpoint label
    # and put tick marks
    nlabs = len(labels)
    step = len(bands.kpoints) / (nlabs - 1)
    for i, lab in enumerate(labels):
        ax1.vlines(i * step, emin, emax, "k")
    ax1.set_xticks([i * step for i in range(nlabs)])
    ax1.set_xticklabels(labels)

    # --------------------------------------------------------------
    # Density of states
    # --------------------------------------------------------------
    # on the x -axis we have the density of states
    # to plot both spin up and down density in the same plot
    # set the xlimit from a negative to positive value
    # ax2.set_title("Density of States")    
    ax2.set_xlim(-10, 10)
    ax2.set_xticklabels([])
    ax2.set_xlabel("Density of States", labelpad=28)
    ax2.set_yticklabels([])
    ax2.grid()
    # line through fermi level, E=0
    ax2.hlines(y=0, xmin=0, xmax=8, color="k", lw=2)

    # plot spd contributions, SPIN UP, positive direction
    ax2.plot(spd_dos["S"].densities[Spin.up],
             run.tdos.energies - run.efermi,
             "r-", label="3s", lw=2)
    ax2.plot(spd_dos["P"].densities[Spin.up],
             run.tdos.energies - run.efermi,
             "g-", label="3p", lw=2)
    ax2.plot(spd_dos["D"].densities[Spin.up],
             run.tdos.energies - run.efermi,
             "b-", label="3d", lw=2)

    # plot spd contribution, SPIN DOWN, negative direction
    ax2.plot(-spd_dos["S"].densities[Spin.down],
             run.tdos.energies - run.efermi,
             "r-", label="3s", lw=2)
    ax2.plot(-spd_dos["P"].densities[Spin.down],
             run.tdos.energies - run.efermi,
             "g-", label="3p", lw=2)
    ax2.plot(-spd_dos["D"].densities[Spin.down],
             run.tdos.energies - run.efermi,
             "b-", label="3d", lw=2)

    # plot total dos
    total_dos = run.tdos.densities[Spin.up] + run.tdos.densities[Spin.down]
    ax2.fill_between(total_dos, 0, run.tdos.energies - run.efermi,
                     color=(0.7, 0.7, 0.7), facecolor=(0.7, 0.7, 0.7))
    ax2.plot(total_dos, run.tdos.energies - run.efermi,
             color=(0.3, 0.3, 0.3), label="total DOS")
    # set legend
    # ax2.legend(fancybox=True, shadow=True, prop={'size':18})
    plt.subplots_adjust(wspace=0)

    # save plot as pdf
    # plt.show()
    plt.savefig("band_dos_up_down.pdf", format="pdf")
