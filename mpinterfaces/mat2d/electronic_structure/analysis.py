from __future__ import print_function, division, unicode_literals

import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun, Locpot, VolumetricData
from pymatgen.io.vasp.inputs import Incar
from pymatgen.electronic_structure.plotter import BSPlotter, BSPlotterProjected
from pymatgen.electronic_structure.core import Spin

from mpinterfaces.utils import is_converged

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


def get_band_edges():
    """
    Calculate the band edge locations relative to the vacuum level
    for a semiconductor.

    Returns:
        edges (dict): {'up_cbm': , 'up_vbm': , 'dn_cbm': , 'dn_vbm': , 'efermi'}
    """
    # Vacuum level energy from LOCPOT.
    locpot = Locpot.from_file('LOCPOT')
    evac = max(locpot.get_average_along_axis(2))

    vasprun = Vasprun('vasprun.xml')
    bs = vasprun.get_band_structure()
    eigenvals = vasprun.eigenvalues
    efermi = vasprun.efermi - evac

    if bs.is_spin_polarized:
        print(eigenvals[Spin.up])
        print([e[0]-evac for e in eigenvals[Spin.up][0]])
        up_cbm = min(
            [min([e[0] for e in eigenvals[Spin.up][i] if not e[1]])
             for i in range(len(eigenvals[Spin.up]))]) - evac
        up_vbm = max(
            [max([e[0] for e in eigenvals[Spin.up][i] if e[1]])
             for i in range(len(eigenvals[Spin.up]))]) - evac
        dn_cbm = min(
            [min([e[0] for e in eigenvals[Spin.down][i] if not e[1]])
             for i in range(len(eigenvals[Spin.down]))]) - evac
        dn_vbm = max(
            [max([e[0] for e in eigenvals[Spin.down][i] if e[1]])
             for i in range(len(eigenvals[Spin.down]))]) - evac
        edges = {'up_cbm': up_cbm, 'up_vbm': up_vbm, 'dn_cbm': dn_cbm,
                 'dn_vbm': dn_vbm, 'efermi': efermi}

    else:
        cbm = bs.get_cbm()['energy'] - evac
        vbm = bs.get_vbm()['energy'] - evac
        edges = {'up_cbm': cbm, 'up_vbm': vbm, 'dn_cbm': cbm, 'dn_vbm': vbm,
                 'efermi': efermi}

    return edges


def plot_band_alignments(directories, run_type='PBE', fmt='pdf'):
    """
    Plot CBM's and VBM's of all compounds together, relative to the band
    edges of H2O.

    Args:
        directories (list): list of the directory paths for materials
            to include in the plot.
        run_type (str): 'PBE' or 'HSE', so that the function knows which
            subdirectory to go into (pbe_bands or hse_bands).
        fmt (str): matplotlib format style. Check the matplotlib
            docs for options.
    """

    if run_type == 'HSE':
        subdirectory = 'hse_bands'
    else:
        subdirectory = 'pbe_bands'

    band_gaps = {}

    for directory in directories:
        sub_dir = os.path.join(directory, subdirectory)
        if is_converged(sub_dir):
            os.chdir(sub_dir)
            band_structure = Vasprun('vasprun.xml').get_band_structure()
            band_gap = band_structure.get_band_gap()

            # Vacuum level energy from LOCPOT.
            locpot = Locpot.from_file('LOCPOT')
            evac = max(locpot.get_average_along_axis(2))

            if not band_structure.is_metal():
                is_direct = band_gap['direct']
                cbm = band_structure.get_cbm()
                vbm = band_structure.get_vbm()

            else:
                cbm = None
                vbm = None
                is_direct = False

            band_gaps[directory] = {'CBM': cbm, 'VBM': vbm,
                                    'Direct': is_direct,
                                    'Metal': band_structure.is_metal(),
                                    'E_vac': evac}

            os.chdir('../../')

    ax = plt.figure(figsize=(16, 10)).gca()

    x_max = len(band_gaps) * 1.315
    ax.set_xlim(0, x_max)

    # Rectangle representing band edges of water.
    ax.add_patch(plt.Rectangle((0, -5.67), height=1.23, width=len(band_gaps),
                 facecolor='#00cc99', linewidth=0))

    ax.text(len(band_gaps) * 1.01, -4.44, r'$\mathrm{H+/H_2}$', size=20,
            verticalalignment='center')
    ax.text(len(band_gaps) * 1.01, -5.67, r'$\mathrm{O_2/H_2O}$', size=20,
            verticalalignment='center')

    x_ticklabels = []

    y_min = -8

    i = 0

    # Nothing but lies.
    are_directs, are_indirects, are_metals = False, False, False

    for compound in [cpd for cpd in directories if cpd in band_gaps]:
        x_ticklabels.append(compound)

        # Plot all energies relative to their vacuum level.
        evac = band_gaps[compound]['E_vac']
        if band_gaps[compound]['Metal']:
            cbm = -8
            vbm = -2
        else:
            cbm = band_gaps[compound]['CBM']['energy'] - evac
            vbm = band_gaps[compound]['VBM']['energy'] - evac

        # Add a box around direct gap compounds to distinguish them.
        if band_gaps[compound]['Direct']:
            are_directs = True
            linewidth = 5
        elif not band_gaps[compound]['Metal']:
            are_indirects = True
            linewidth = 0

        # Metals are grey.
        if band_gaps[compound]['Metal']:
            are_metals = True
            linewidth = 0
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

    ax.set_ylim(y_min, 0)

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

    # Who needs axes?
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    ax.set_ylabel('eV', family='serif', size=24)

    if fmt == "None":
        return ax
    else:
        plt.savefig('band_alignments.{}'.format(fmt), transparent=True)
    plt.close()


def plot_local_potential(axis=2, ylim=(-20, 0), fmt='pdf'):
    """
    Plot data from the LOCPOT file along any of the 3 primary axes.
    Useful for determining surface dipole moments and electric
    potentials on the interior of the material.

    Args:
        axis (int): 0 = x, 1 = y, 2 = z
        ylim (tuple): minimum and maximum potentials for the plot's y-axis.
        fmt (str): matplotlib format style. Check the matplotlib docs
            for options.
    """

    ax = plt.figure(figsize=(16, 10)).gca()

    locpot = Locpot.from_file('LOCPOT')
    structure = Structure.from_file('CONTCAR')
    vd = VolumetricData(structure, locpot.data)
    abs_potentials = vd.get_average_along_axis(axis)
    vacuum_level = max(abs_potentials)

    vasprun = Vasprun('vasprun.xml')
    bs = vasprun.get_band_structure()
    if not bs.is_metal():
        cbm = bs.get_cbm()['energy'] - vacuum_level
        vbm = bs.get_vbm()['energy'] - vacuum_level

    potentials = [potential - vacuum_level for potential in abs_potentials]
    axis_length = structure.lattice._lengths[axis]
    positions = np.arange(0, axis_length, axis_length / len(potentials))

    ax.plot(positions, potentials, linewidth=2, color='k')

    ax.set_xlim(0, axis_length)
    ax.set_ylim(ylim[0], ylim[1])

    ax.set_xticklabels(
        [r'$\mathrm{%s}$' % tick for tick in ax.get_xticks()], size=20)
    ax.set_yticklabels(
        [r'$\mathrm{%s}$' % tick for tick in ax.get_yticks()], size=20)
    ax.set_xlabel(r'$\mathrm{\AA}$', size=24)
    ax.set_ylabel(r'$\mathrm{V\/(eV)}$', size=24)

    if not bs.is_metal():
        ax.text(ax.get_xlim()[1], cbm, r'$\mathrm{CBM}$',
                horizontalalignment='right', verticalalignment='bottom',
                size=20)
        ax.text(ax.get_xlim()[1], vbm, r'$\mathrm{VBM}$',
                horizontalalignment='right', verticalalignment='top', size=20)
        ax.fill_between(ax.get_xlim(), cbm, ax.get_ylim()[1],
                        facecolor=plt.cm.jet(0.3), zorder=0, linewidth=0)
        ax.fill_between(ax.get_xlim(), ax.get_ylim()[0], vbm,
                        facecolor=plt.cm.jet(0.7), zorder=0, linewidth=0)

    if fmt == "None":
        return ax
    else:
        plt.savefig('locpot.{}'.format(fmt))
    plt.close()


def plot_band_structure(ylim=(-5, 5), draw_fermi=False, fmt='pdf'):
    """
    Plot a standard band structure with no projections.

    Args:
        ylim (tuple): minimum and maximum potentials for the plot's y-axis.
        draw_fermi (bool): whether or not to draw a dashed line at E_F.
        fmt (str): matplotlib format style. Check the matplotlib docs
            for options.
    """

    vasprun = Vasprun('vasprun.xml')
    efermi = vasprun.efermi
    bsp = BSPlotter(vasprun.get_band_structure('KPOINTS', line_mode=True,
                                               efermi=efermi))
    if fmt == "None":
        return bsp.bs_plot_data()
    else:
        plot = bsp.get_plot(ylim=ylim)
        fig = plot.gcf()
        ax = fig.gca()
        ax.set_xticklabels([r'$\mathrm{%s}$' % t for t in ax.get_xticklabels()])
        ax.set_yticklabels([r'$\mathrm{%s}$' % t for t in ax.get_yticklabels()])
        if draw_fermi:
            ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [0, 0], 'k--')
        plt.savefig('band_structure.{}'.format(fmt), transparent=True)

    plt.close()


def plot_color_projected_bands(ylim=(-5, 5), fmt='pdf'):
    """
    Plot a single band structure where the color of the band indicates
    the elemental character of the eigenvalue.

    Args:
        ylim (tuple): minimum and maximum energies for the plot's
            y-axis.
        fmt (str): matplotlib format style. Check the matplotlib
            docs for options.
    """

    vasprun = Vasprun('vasprun.xml', parse_projected_eigen=True)
    bs = vasprun.get_band_structure('KPOINTS', line_mode=True)
    bspp = BSPlotterProjected(bs)
    ax = bspp.get_elt_projected_plots_color().gcf().gca()
    ax.set_xticklabels([r'$\mathrm{%s}$' % t for t in ax.get_xticklabels()])
    ax.set_yticklabels([r'$\mathrm{%s}$' % t for t in ax.get_yticklabels()])
    ax.set_ylim(ylim)
    if fmt == "None":
        return ax
    else:
        plt.savefig('color_projected_bands.{}'.format(fmt))
    plt.close()


def plot_elt_projected_bands(ylim=(-5, 5), fmt='pdf'):
    """
    Plot separate band structures for each element where the size of the
    markers indicates the elemental character of the eigenvalue.

    Args:
        ylim (tuple): minimum and maximum energies for the plot's
            y-axis.
        fmt (str): matplotlib format style. Check the matplotlib
            docs for options.
    """

    vasprun = Vasprun('vasprun.xml', parse_projected_eigen=True)
    bs = vasprun.get_band_structure('KPOINTS', line_mode=True)
    bspp = BSPlotterProjected(bs)
    ax = bspp.get_elt_projected_plots(ylim=ylim).gcf().gca()
    ax.set_xticklabels([r'$\mathrm{%s}$' % t for t in ax.get_xticklabels()])
    ax.set_yticklabels([r'$\mathrm{%s}$' % t for t in ax.get_yticklabels()])
    if fmt == "None":
        return ax
    else:
        plt.savefig('elt_projected_bands.{}'.format(fmt))
    plt.close()


def plot_orb_projected_bands(orbitals, fmt='pdf', ylim=(-5, 5)):
    """
    Plot a separate band structure for each orbital of each element in
    orbitals.

    Args:
        orbitals (dict): dictionary of the form
            {element: [orbitals]},
            e.g. {'Mo': ['s', 'p', 'd'], 'S': ['p']}
        ylim (tuple): minimum and maximum energies for the plot's
            y-axis.
        fmt (str): matplotlib format style. Check the matplotlib
            docs for options.
    """

    vasprun = Vasprun('vasprun.xml', parse_projected_eigen=True)
    bs = vasprun.get_band_structure('KPOINTS', line_mode=True)
    bspp = BSPlotterProjected(bs)
    ax = bspp.get_projected_plots_dots(orbitals, ylim=ylim).gcf().gca()
    ax.set_xticklabels([r'$\mathrm{%s}$' % t for t in ax.get_xticklabels()])
    ax.set_yticklabels([r'$\mathrm{%s}$' % t for t in ax.get_yticklabels()])
    if fmt == "None":
        return ax
    else:
        plt.savefig('orb_projected_bands.{}'.format(fmt))
    plt.close()


def get_effective_mass():
    """
    This function is in a beta stage, and its results are not
    guaranteed to be useful.

    Finds effective masses from a band structure, using parabolic
    fitting to determine the band curvature at the CBM
    for electrons and at the VBM for holes. This curvature enters
    the equation m* = (hbar)**2 / (d^2E/dk^2).

    To consider anisotropy, the k-space directions to the left and right
    of the CBM/VBM in the band diagram are returned separately.

    *NOTE* Only works for semiconductors and linemode calculations (non-
           spin polarized).
           >30 k-points per string recommended to obtain
           reliable curvatures.

    *NOTE* The parabolic fit can be quite sensitive to the number of
           k-points fit to, so it might be worthwhile adjusting N_KPTS
           to obtain some sense of the error bar.

    TODO: Warn user if CBM/VBM is at the edge of the diagram, and which
          direction (either left or right) was not actually fit to.
          Until fixed, this (most likely) explains any negative masses
          returned.

    Returns:
        Dictionary of the form
            {'electron': {'left': e_m_eff_l, 'right': e_m_eff_r},
             'hole': {'left': h_m_eff_l, 'right': h_m_eff_r}}
            where 'left' and 'right' indicate the reciprocal
            directions to the left and right of the extremum in the
            band structure.
    """

    H_BAR = 6.582119514e-16  # eV*s
    M_0 = 9.10938356e-31  # kg
    N_KPTS = 6  # Number of k-points included in the parabola.

    spin_up = Spin(1)

    band_structure = Vasprun('vasprun.xml').get_band_structure()

    # Locations of CBM and VBM in band_structure.bands
    cbm_band_index = band_structure.get_cbm()['band_index'][spin_up][0]
    cbm_kpoint_index = band_structure.get_cbm()['kpoint_index'][0]

    vbm_band_index = band_structure.get_vbm()['band_index'][spin_up][0]
    vbm_kpoint_index = band_structure.get_vbm()['kpoint_index'][0]

    k = {'electron': {'left': [], 'right': []},
         'hole': {'left': [], 'right': []}}
    E = {'electron': {'left': [], 'right': []},
         'hole': {'left': [], 'right': []}}

    e_ref_coords = band_structure.kpoints[cbm_kpoint_index]._ccoords
    h_ref_coords = band_structure.kpoints[vbm_kpoint_index]._ccoords

    for n in range(-N_KPTS, 1):
        e_coords = band_structure.kpoints[cbm_kpoint_index + n]._ccoords
        h_coords = band_structure.kpoints[vbm_kpoint_index + n]._ccoords

        k['electron']['left'].append(
            ((e_coords[0] - e_ref_coords[0])**2 +
             (e_coords[1] - e_ref_coords[1])**2 +
             (e_coords[2] - e_ref_coords[2])**2)**0.5
            )
        k['hole']['left'].append(
            ((h_coords[0] - h_ref_coords[0])**2 +
             (h_coords[1] - h_ref_coords[1])**2 +
             (h_coords[2] - h_ref_coords[2])**2)**0.5
            )

        e_energy = band_structure.bands[
            spin_up][cbm_band_index][cbm_kpoint_index + n]
        h_energy = band_structure.bands[
            spin_up][vbm_band_index][vbm_kpoint_index + n]

        E['electron']['left'].append(e_energy)
        E['hole']['left'].append(h_energy)

    for n in range(1, 1 + N_KPTS):
        e_coords = band_structure.kpoints[cbm_kpoint_index + n]._ccoords
        h_coords = band_structure.kpoints[vbm_kpoint_index + n]._ccoords

        k['electron']['right'].append(
            ((e_coords[0] - e_ref_coords[0])**2 +
             (e_coords[1] - e_ref_coords[1])**2 +
             (e_coords[2] - e_ref_coords[2])**2)**0.5
            )
        k['hole']['right'].append(
            ((h_coords[0] - h_ref_coords[0])**2 +
             (h_coords[1] - h_ref_coords[1])**2 +
             (h_coords[2] - h_ref_coords[2])**2)**0.5
            )

        e_energy = band_structure.bands[
            spin_up][cbm_band_index][cbm_kpoint_index + n]
        h_energy = band_structure.bands[
            spin_up][vbm_band_index][vbm_kpoint_index + n]

        E['electron']['right'].append(e_energy)
        E['hole']['right'].append(h_energy)

    # 2nd order fits
    e_l_fit = np.poly1d(
        np.polyfit(k['electron']['left'], E['electron']['left'], 2))
    e_r_fit = np.poly1d(
        np.polyfit(k['electron']['right'], E['electron']['right'], 2))
    h_l_fit = np.poly1d(
        np.polyfit(k['hole']['left'], E['hole']['left'], 2))
    h_r_fit = np.poly1d(
        np.polyfit(k['hole']['right'], E['hole']['right'], 2))

    # Curvatures
    e_l_curvature = e_l_fit.deriv().deriv()[0]
    e_r_curvature = e_r_fit.deriv().deriv()[0]
    h_l_curvature = h_l_fit.deriv().deriv()[0]
    h_r_curvature = h_r_fit.deriv().deriv()[0]

    # Unit conversion
    e_m_eff_l = 10 * ((H_BAR ** 2) / e_l_curvature) / M_0
    e_m_eff_r = 10 * ((H_BAR ** 2) / e_r_curvature) / M_0
    h_m_eff_l = -10 * ((H_BAR ** 2) / h_l_curvature) / M_0
    h_m_eff_r = -10 * ((H_BAR ** 2) / h_r_curvature) / M_0

    return {'electron': {'left': e_m_eff_l, 'right': e_m_eff_r},
            'hole': {'left': h_m_eff_l, 'right': h_m_eff_r}}


def plot_density_of_states(xlim=(-10, 5), ylim=(-1.5, 1.5), fmt='pdf'):
    """
    Plots the density of states from the DOSCAR in the cwd. Plots
    spin up in red, down in green, and the sum in black. Efermi = 0.

    Args:
        xlim (tuple): minimum and maximum energies for the plot's
            x-axis.
        ylim (tuple): minimum and maximum for the plot's
            y-axis.
        fmt (str): matplotlib format style. Check the matplotlib
            docs for options.
    """

    efermi = Vasprun('vasprun.xml').efermi
    dos_lines = open ('DOSCAR').readlines()

    x, up, down = [], [], []
    nedos = Incar.from_file('INCAR').as_dict()['NEDOS'] - 1

    for line in dos_lines[6:6+nedos]:
            split_line = line.split()
            x.append(float(split_line[0]) - efermi)
            up.append(float(split_line[1]))
            down.append(-float(split_line[2]))

    x, up, down = np.array(x), np.array(up), np.array(down)
    sum = up + down

    ax = plt.figure().gca()
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])

    ax.set_xlabel(r'$\mathrm{E\/(eV)}$')
    ax.set_ylabel(r'$\mathrm{Density\/of\/States$')
    ax.set_xticklabels([r'$\mathrm{%s}$' % t for t in ax.get_xticklabels()])
    ax.set_yticklabels([r'$\mathrm{%s}$' % t for t in ax.get_yticklabels()])

    ax.plot(x, up, color='red' )
    ax.plot(x, down, color='green')
    ax.plot(x, sum, color='black' )
    if fmt is not None:
        plt.savefig('density_of_states.{}'.format(fmt))
    else:
        return ax

    plt.close()


def get_fermi_velocities():
    """
    Calculates the fermi velocity of each band that crosses the fermi
    level, according to v_F = dE/(h_bar*dk).

    Returns:
        fermi_velocities (list). The absolute values of the
            adjusted slopes of each band, in Angstroms/s.
    """

    vr = Vasprun('vasprun.xml')
    # eigenvalues = vr.eigenvalues
    bs = vr.get_band_structure()
    bands = bs.bands
    kpoints = bs.kpoints
    efermi = bs.efermi
    h_bar = 6.582e-16  # eV*s

    fermi_bands = []
    for spin in bands:
        for i in range(len(bands[spin])):
            if max(bands[spin][i]) > efermi > min(bands[spin][i]):
                fermi_bands.append(bands[spin][i])

    fermi_velocities = []
    for band in fermi_bands:
        for i in range(len(band)-1):
            if (band[i] < efermi < band[i+1]) or (band[i] > efermi > band[i+1]):
                dk = np.sqrt((kpoints[i+1].cart_coords[0]
                              - kpoints[i].cart_coords[0])**2
                             + (kpoints[i+1].cart_coords[1]
                                - kpoints[i].cart_coords[1])**2)
                v_f = abs((band[i+1] - band[i]) / (h_bar * dk))
                fermi_velocities.append(v_f)

    return fermi_velocities  # Values are in Angst./s


def find_dirac_nodes():
    """
    Look for band crossings near (within `tol` eV) the Fermi level.

    Returns:
        boolean. Whether or not a band crossing occurs at or near
            the fermi level.
    """

    vasprun = Vasprun('vasprun.xml')
    dirac = False
    if vasprun.get_band_structure().get_band_gap()['energy'] < 0.1:
        efermi = vasprun.efermi
        bsp = BSPlotter(vasprun.get_band_structure('KPOINTS', line_mode=True,
                                                   efermi=efermi))
        bands = []
        data = bsp.bs_plot_data(zero_to_efermi=True)
        for d in range(len(data['distances'])):
            for i in range(bsp._nb_bands):
                x = data['distances'][d],
                y = [data['energy'][d][str(Spin.up)][i][j]
                     for j in range(len(data['distances'][d]))]
                band = [x, y]
                bands.append(band)

        considered = []
        for i in range(len(bands)):
            for j in range(len(bands)):
                if i != j and (j, i) not in considered:
                    considered.append((j, i))
                    for k in range(len(bands[i][0])):
                        if ((-0.1 < bands[i][1][k] < 0.1) and
                                (-0.1 < bands[i][1][k] - bands[j][1][k] < 0.1)):
                            dirac = True
    return dirac


def plot_spin_texture(inner_index, outer_index, center=(0, 0), fmt='pdf'):
    """
    Create six plots- one for the spin texture in x, y, and z in
    each of two bands: an inner band and an outer band. For
    Rashba spin-splitting, these two bands should be the two that
    have split.

    Args:
        inner_index (int): indices of the two spin-split bands.
        outer_index (int): indices of the two spin-split bands.
        center (tuple): coordinates of the center of the splitting
            (where the bands cross). Defaults to Gamma.
        fmt: matplotlib format style. Check the matplotlib
            docs for options.
    """

    procar_lines = open("PROCAR").readlines()

    data = procar_lines[1].split()
    n_kpts = int(data[3])
    n_bands = int(data[7])
    n_ions = int(data[11])

    # These numbers, along with almost everything else in this
    # function, are magical. Don't touch them.
    band_step = (n_ions + 1) * 4 + 4
    k_step = n_bands * band_step + 3

    kpoints = []
    spin_textures = {'inner': {'x': [], 'y': [], 'z': []},
                     'outer': {'x': [], 'y': [], 'z': []}}

    for n in range(n_kpts):
        for var in ['x', 'y', 'z']:
            spin_textures['inner'][var].append(0)
            spin_textures['outer'][var].append(0)

    i = 3
    j = 0
    while i < len(procar_lines):
        kpoints.append([float(procar_lines[i][18:29]) - center[0],
                        float(procar_lines[i][29:40]) - center[1]])
        spin_textures['inner']['x'][j] += float(
            procar_lines[i+(4+(n_ions+1)*2)+inner_index*band_step].split()[-1])
        spin_textures['inner']['y'][j] += float(
            procar_lines[i+(4+(n_ions+1)*3)+inner_index*band_step].split()[-1])
        spin_textures['inner']['z'][j] += float(
            procar_lines[i+(4+(n_ions+1)*4)+inner_index*band_step].split()[-1])
        spin_textures['outer']['x'][j] += float(
            procar_lines[i+(4+(n_ions+1)*2)+outer_index*band_step].split()[-1])
        spin_textures['outer']['y'][j] += float(
            procar_lines[i+(4+(n_ions+1)*3)+outer_index*band_step].split()[-1])
        spin_textures['outer']['z'][j] += float(
            procar_lines[i+(4+(n_ions+1)*4)+outer_index*band_step].split()[-1])
        i += k_step
        j += 1

    for branch in spin_textures:
        for vector in spin_textures[branch]:
            print('plotting {}_{}.{}'.format(branch, vector, fmt))
            ax = plt.subplot(111, projection='polar')

            raw = [spin_textures[branch][vector][k] for k in range(len(kpoints))]
            minimum = min(raw)
            maximum = max(raw) - minimum

            r_max = max([np.sqrt(kpt[0]**2 + kpt[1]**2) for kpt in kpoints])

            for l in range(len(kpoints)):
                if kpoints[l][0] == 0 and kpoints[l][1] > 0:
                    theta = np.pi / 2.0
                elif kpoints[l][0] == 0:
                    theta = 3.0 * np.pi / 2.0
                elif kpoints[l][0] < 0:
                    theta = np.pi + np.arctan(kpoints[l][1] / kpoints[l][0])
                else:
                    theta = np.arctan(kpoints[l][1] / kpoints[l][0])
                r = np.sqrt(kpoints[l][0]**2 + kpoints[l][1]**2)
                if r == 0:
                    w = 0
                else:
                    w = r_max*0.07/r
                ax.add_patch(
                    plt.Rectangle(
                        (theta, r), width=w, height=r_max*0.07,
                        color=plt.cm.rainbow(
                            (spin_textures[branch][vector][l]-minimum)/maximum
                        )
                    )
                )
            ax.plot(0, 0, linewidth=0, marker='o', color='k', markersize=18)
            ax.set_rmax(r_max)
            plt.axis('off')
            plt.savefig('{}_{}.{}'.format(branch, vector, fmt))
            plt.close()
