import os

from pymatgen.io.vasp.outputs import Vasprun

import matplotlib.pyplot as plt


def get_band_structures(directories):

    band_gaps = {}
    for directory in directories:
        os.chdir(directory)
        vasprun = Vasprun('vasprun.xml')
        band_gap = vasprun.get_band_structure().get_band_gap()

        transition = band_gap['transition'].split('-')

        if transition[0] == transition[1]:
            is_direct = True
        else:
            is_direct = False

        if band_gap['energy']:
            cbm = vasprun.get_band_structure().get_cbm()
            vbm = vasprun.get_band_structure().get_vbm()
            band_gaps[directory] = {'CBM': cbm, 'VBM': vbm,
                                    'Direct': is_direct}
        else:
            band_gaps[directory] = False

    return band_gaps


def plot_band_alignments(band_gaps):

    ax = plt.figure(figsize=(10, 10)).gca()
    ax.set_xlim(len(band_gaps))

    x_ticklabels = []
    vbms = []
    i = 0
    for compound in band_gaps:
        x_ticklabels.append(compound)
        cbm = band_gaps[compound]['CBM']
        vbm = band_gaps[compound]['VBM']
        vbms.append(vbm)
        if band_gaps[compound]['Direct']:
            linewidth = 3
        else:
            linewidth = 0

        ax.add_patch(plt.Rectangle((i, cbm), height=-cbm, width=0.9,
                                   facecolor="#0066cc", linewidth=linewidth))
        i += 1

    ax.set_ylim(min(vbms)-0.5, 0)

    ax.plot([0, i], [-1.9, -1.9], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -1.94, r'$\mathrm{CO_2+\/e^-\/\rightarrow\/CO^-_2}$',
            size=20)

    ax.plot([0, i], [-0.61, -0.61], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -0.7,
            r'$\mathrm{CO_2\/+\/2H^+\/+\/2e^-\/\rightarrow\/HCO_2H}$', size=20)

    ax.plot([0, i], [-0.53, -0.53], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -0.59,
            r'$\mathrm{CO_2\/+\/2H^+\/+\/2e^-\/\rightarrow\/CO\/+\/H_2O}$',
            size=20)

    ax.plot([0, i], [-0.48, -0.48], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -0.48,
            r'$\mathrm{CO_2\/+\/4H^+\/+\/4e^-\/\rightarrow\/HCHO\/+\/H_2O}$',
            size=20)

    ax.plot([0, i], [-0.38, -0.38], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -0.37,
            r'$\mathrm{CO_2\/+\/6H^+\/+\/6e^-\/\rightarrow\/CH_3OH\/+\/H_2O}$',
            size=20)

    ax.plot([0, 16], [-0.24, -0.24], color='k', alpha=0.6, linewidth=4)
    ax.text(i + 0.5, -0.24,
            r'$\mathrm{CO_2\/+\/8H^+\/+\/8e^-\/\rightarrow\/CH_4\/+\/2H_2O}$',
            size=20)

    ax.set_xticks(range(i))
    ax.set_xticklabels(x_ticklabels, family='serif', size=20, rotation=60)
    ax.set_yticklabels(ax.get_yticks(), family='serif', size=20)

    ax.add_patch(plt.Rectangle((i + 1, min(vbms) + 0.25), width=3, height=0.2,
                               facecolor='w', linewidth=3))
    ax.text(i + 2.5, min(vbms) + 0.35, 'Direct', family='serif', color='k',
            size=20, horizontalalignment='center', verticalalignment='center')
    ax.add_patch(plt.Rectangle((i + 5, min(vbms) + 0.25), width=3, height=0.2,
                               facecolor='w', linewidth=0))
    ax.text(i + 6.5, min(vbms) + 0.35, 'Indirect', family='serif', size=20,
            color='k', horizontalalignment='center',
            verticalalignment='center')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    plt.savefig('band_alignments.pdf', transparent=True)
