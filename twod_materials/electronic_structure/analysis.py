import os

import numpy as np

from pymatgen.io.vasp.outputs import Vasprun, Locpot
from pymatgen.electronic_structure.plotter import BSPlotter

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from twod_materials.utils import is_converged


def plot_band_alignments(directories):
    """
    Plot CBM's and VBM's of all compounds together, relative to the band
    edges of H2O.
    """

    band_gaps = {}

    for directory in directories:
        if is_converged('{}/bandstructure'.format(directory)):
            os.chdir('{}/bandstructure'.format(directory))
            vasprun = Vasprun('vasprun.xml')
            band_gap = vasprun.get_band_structure().get_band_gap()

            # Vacuum level energy from LOCPOT.
            locpot = Locpot.from_file('LOCPOT')
            evac = locpot.get_average_along_axis(2)[-5]

            try:
                transition = band_gap['transition'].split('-')

                if transition[0] == transition[1]:
                    is_direct = True
                else:
                    is_direct = False

                is_metal = False
                cbm = vasprun.get_band_structure().get_cbm()
                vbm = vasprun.get_band_structure().get_vbm()

            except AttributeError:
                cbm = None
                vbm = None
                is_metal = True
                is_direct = False

            band_gaps[directory] = {'CBM': cbm, 'VBM': vbm,
                                    'Direct': is_direct, 'Metal': is_metal,
                                    'E_vac': evac}

            os.chdir('../../')

    ax = plt.figure(figsize=(16, 10)).gca()

    x_max = len(band_gaps)*1.315
    ax.set_xlim(0, x_max)

    # Rectangle representing band edges of water.
    ax.add_patch(plt.Rectangle((0, -5.67), height=1.23, width=len(band_gaps),
                 facecolor='#00cc99', linewidth=0))
    ax.text(len(band_gaps)*1.01, -4.44, r'$\mathrm{H+/H_2}$', size=20,
            verticalalignment='center')
    ax.text(len(band_gaps)*1.01, -5.67, r'$\mathrm{O2/H_2O}$', size=20,
            verticalalignment='center')

    x_ticklabels = []
    vbms = []
    for compound in band_gaps:
        vbms.append(band_gaps[compound]['VBM']['energy'])

    y_min = -8

    i = 0

    # Nothing but lies.
    are_directs, are_indirects, are_metals = False, False, False

    for compound in band_gaps:
        x_ticklabels.append(compound)

        # Plot all energies relative to their vacuum level.
        evac = band_gaps[compound]['E_vac']
        cbm = band_gaps[compound]['CBM']['energy'] - evac
        vbm = band_gaps[compound]['VBM']['energy'] - evac

        # Add a box around direct gap compounds to distinguish them.
        if band_gaps[compound]['Direct']:
            are_directs = True
            linewidth = 5
        else:
            are_indirects = True
            linewidth = 0

        # Metals are grey.
        if band_gaps[compound]['Metal']:
            are_metals = True
            color_code = '#404040'
        else:
            color_code = '#002b80'

        # CBM
        ax.add_patch(plt.Rectangle((i, cbm), height=-cbm, width=0.8,
                                   facecolor=color_code, linewidth=linewidth,
                                   edgecolor="#e68a00"))
        # VBM
        ax.add_patch(plt.Rectangle((i, y_min),
                                   height=(vbm - y_min), width=0.8,
                                   facecolor=color_code, linewidth=linewidth,
                                   edgecolor="#e68a00"))

        i += 1

    ax.set_ylim(y_min, -2)

    # Set tick labels
    ax.set_xticks([n + 0.4 for n in range(i)])
    ax.set_xticklabels(x_ticklabels, family='serif', size=20, rotation=60)
    ax.set_yticklabels(ax.get_yticks(), family='serif', size=20)

    # Add a legend
    height = y_min
    if are_directs:
        ax.add_patch(plt.Rectangle((i*1.165, height), width=i*0.15,
                                   height=(-y_min*0.1), facecolor='#002b80',
                                   edgecolor='#e68a00', linewidth=5))
        ax.text(i*1.24, height - y_min * 0.05, 'Direct', family='serif',
                color='w', size=20, horizontalalignment='center',
                verticalalignment='center')
        height -= y_min * 0.15

    if are_indirects:
        ax.add_patch(plt.Rectangle((i*1.165, height), width=i*0.15,
                                   height=(-y_min*0.1), facecolor='#002b80',
                                   linewidth=0))
        ax.text(i*1.24, height - y_min * 0.05, 'Indirect', family='serif',
                size=20, color='w', horizontalalignment='center',
                verticalalignment='center')
        height -= y_min * 0.15

    if are_metals:
        ax.add_patch(plt.Rectangle((i*1.165, height), width=i*0.15,
                                   height=(-y_min*0.1), facecolor='#404040',
                                   linewidth=0))
        ax.text(i*1.24, height - y_min * 0.05, 'Metal', family='serif',
                size=20, color='w', horizontalalignment='center',
                verticalalignment='center')

    # Axes are hideous.
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    plt.savefig('band_alignments.pdf', transparent=True)


def plot_band_structure():
    vasprun = Vasprun('vasprun.xml')
    bsp = BSPlotter(vasprun.get_band_structure('KPOINTS', line_mode=True,
                                               efermi=vasprun.efermi))
    bsp.save_plot('band_structure.eps', ylim=(-5, 5))


def get_effective_mass():
    """
    Returns qualitative, but not exactly quantitatively accurate
    effective masses from a band structure. Uses parabolic
    fitting to determine the band curvature at the CBM
    for electrons and at the VBM for holes. This curvature enters
    the equation m* = (hbar)**2 / (d^2E/dk^2).

    *NOTE* Only works for semiconductors and linemode calculations.
           >30 k-points per string recommended to obtain
           reliable curvatures.
    """

    H_BAR = 6.582119514e-16  # eV*s
    M_0 = 9.10938356e-31  # kg

    band_structure = Vasprun('vasprun.xml').get_band_structure()

    cbm_band_index = band_structure.get_cbm()['band_index'][1][0]
    cbm_kpoint_index = band_structure.get_cbm()['kpoint_index'][0]

    vbm_band_index = band_structure.get_vbm()['band_index'][1][0]
    vbm_kpoint_index = band_structure.get_vbm()['kpoint_index'][0]

    e_k = {'left': [], 'right': []}
    e_E = {'left': [], 'right': []}
    h_k = {'left': [], 'right': []}
    h_E = {'left': [], 'right': []}

    e_ref_coords = band_structure.kpoints[cbm_kpoint_index]._ccoords
    h_ref_coords = band_structure.kpoints[vbm_kpoint_index]._ccoords

    for n in range(-6, 1):
        e_coords = band_structure.kpoints[cbm_kpoint_index + n]._ccoords
        h_coords = band_structure.kpoints[vbm_kpoint_index + n]._ccoords

        e_k['left'].append(
            ((e_coords[0] - e_ref_coords[0])**2 +
             (e_coords[1] - e_ref_coords[1])**2 +
             (e_coords[2] - e_ref_coords[2])**2)**0.5
            )
        h_k['left'].append(
            ((h_coords[0] - h_ref_coords[0])**2 +
             (h_coords[1] - h_ref_coords[1])**2 +
             (h_coords[2] - h_ref_coords[2])**2)**0.5
            )

        e_energy = band_structure.bands[1][cbm_band_index][cbm_kpoint_index + n]
        h_energy = band_structure.bands[1][vbm_band_index][vbm_kpoint_index + n]

        e_E['left'].append(e_energy)
        h_E['left'].append(h_energy)

    for n in range(1, 7):
        e_coords = band_structure.kpoints[cbm_kpoint_index + n]._ccoords
        h_coords = band_structure.kpoints[vbm_kpoint_index + n]._ccoords

        e_k['right'].append(
            ((e_coords[0] - e_ref_coords[0])**2 +
             (e_coords[1] - e_ref_coords[1])**2 +
             (e_coords[2] - e_ref_coords[2])**2)**0.5
            )
        h_k['right'].append(
            ((h_coords[0] - h_ref_coords[0])**2 +
             (h_coords[1] - h_ref_coords[1])**2 +
             (h_coords[2] - h_ref_coords[2])**2)**0.5
            )

        e_energy = band_structure.bands[1][cbm_band_index][cbm_kpoint_index + n]
        h_energy = band_structure.bands[1][vbm_band_index][vbm_kpoint_index + n]

        e_E['right'].append(e_energy)
        h_E['right'].append(h_energy)

    e_l_fit = np.polyfit(e_k['left'], e_E['left'], 2)
    e_r_fit = np.polyfit(e_k['right'], e_E['right'], 2)
    h_l_fit = np.polyfit(h_k['left'], h_E['left'], 2)
    h_r_fit = np.polyfit(h_k['right'], h_E['right'], 2)

    e_l_function = np.poly1d(e_l_fit)
    e_r_function = np.poly1d(e_r_fit)
    h_l_function = np.poly1d(h_l_fit)
    h_r_function = np.poly1d(h_r_fit)

    e_l_curvature = e_l_function.deriv().deriv()[0]
    e_r_curvature = e_r_function.deriv().deriv()[0]
    h_l_curvature = h_l_function.deriv().deriv()[0]
    h_r_curvature = h_r_function.deriv().deriv()[0]

    e_m_eff_l = 10 * ((H_BAR ** 2) / e_l_curvature) / M_0
    e_m_eff_r = 10 * ((H_BAR ** 2) / e_r_curvature) / M_0
    h_m_eff_l = -10 * ((H_BAR ** 2) / h_l_curvature) / M_0
    h_m_eff_r = -10 * ((H_BAR ** 2) / h_r_curvature) / M_0

    return {'electron': [e_m_eff_l, e_m_eff_r], 'hole': [h_m_eff_l, h_m_eff_r]}

