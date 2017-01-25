from __future__ import print_function, division, unicode_literals

import os

import operator

from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.matproj.rest import MPRester

from monty.serialization import loadfn

from twod_materials.utils import is_converged

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


try:
    MPR = MPRester(
        loadfn(os.path.join(os.path.expanduser('~'), 'config.yaml'))['mp_api']
        )
except IOError:
    try:
        MPR = MPRester(
            os.environ['MP_API']
            )
    except KeyError:
        raise ValueError('No Materials Project API key found. Please check'
                         ' that your ~/config.yaml contains the field'
                         ' mp_api: your_api_key')


def get_competing_phases(directories):
    """
    Collect the species to which the 2D materials might decompose to.
    Since a lot of 2D materials with similar compositions will have the
    same competing phases, duplicates aren't counted.
    """

    total_competing_phases = []

    for directory in directories:
        os.chdir(directory)
        composition = Structure.from_file('POSCAR').composition
        try:
            energy = Vasprun('vasprun.xml').final_energy
        except:
            energy = 100
        my_entry = ComputedEntry(composition, energy)  # 2D material
        entries = MPR.get_entries_in_chemsys(
            [elt.symbol for elt in composition]
            )

        entries.append(my_entry)  # 2D material

        pda = PDAnalyzer(PhaseDiagram(entries))
        decomp = pda.get_decomp_and_e_above_hull(my_entry, allow_negative=True)
        competing_phases = [
            (entry.composition.reduced_formula,
             entry.entry_id) for entry in decomp[0]
            ]

        # Keep a running list of all unique competing phases, since in
        # high throughput 2D searches there is usually some overlap in
        # competing phases for different materials.
        for specie in competing_phases:
            if specie not in total_competing_phases:
                total_competing_phases.append(specie)
        os.chdir('../')

    return total_competing_phases


def get_hull_distances(directories):

    hull_distances = {}
    finished_competitors = {}

    # Determine which competing phases have been relaxed in the current
    # framework and store them in a dictionary ({formula: entry}).
    if os.path.isdir('all_competitors'):
        os.chdir('all_competitors')
        for comp_dir in [
            dir for dir in os.listdir(os.getcwd()) if os.path.isdir(dir) and
            is_converged(dir)
                ]:
            vasprun = Vasprun('{}/vasprun.xml'.format(comp_dir))
            composition = vasprun.final_structure.composition
            energy = vasprun.final_energy
            finished_competitors[comp_dir] = ComputedEntry(composition, energy)
        os.chdir('../')

    for directory in directories:
        os.chdir(directory)
        composition = Structure.from_file('POSCAR').composition
        try:
            energy = Vasprun('vasprun.xml').final_energy
        except:
            energy = 100
        my_entry = ComputedEntry(composition, energy)  # 2D material
        entries = MPR.get_entries_in_chemsys(
            [elt.symbol for elt in composition]
            )

        # If the energies of competing phases have been calculated in
        # the current framework, put them in the phase diagram instead
        # of the MP energies.
        for i in range(len(entries)):
            formula = entries[i].composition.reduced_formula
            if formula in finished_competitors:
                entries[i] = finished_competitors[formula]
            else:
                entries[i] = ComputedEntry(entries[i].composition, 100)

        entries.append(my_entry)  # 2D material

        pda = PDAnalyzer(PhaseDiagram(entries))
        decomp = pda.get_decomp_and_e_above_hull(my_entry, allow_negative=True)

        hull_distances[composition.reduced_formula] = decomp[1]
        os.chdir('../')

    return hull_distances


def plot_hull_distances(hull_distances, fmt='pdf'):
    """
    Create a bar graph of the formation energies of the 2D materials.
    """

    hsize = 12 + (len(hull_distances) - 4) / 3
    ax = plt.figure(figsize=(hsize, 10)).gca()
    ax.set_ylim(0, 700)
    ax.set_xlim(0, len(hull_distances))

    x_ticklabels = []
    i = 0
    for compound in sorted(
            hull_distances.items(), key=operator.itemgetter(1)):

        proper_formula = ''
        for char in compound[0]:
            try:
                int(char)
                proper_formula += '_{}'.format(char)
            except ValueError:
                proper_formula += char

        x_ticklabels.append(r'$\mathrm{%s}$' % proper_formula)
        hull_distance = hull_distances[compound[0]] * 1000

        # Good chance of stability
        if hull_distance < 100:
            color_code = 0.5

        # Decent chance of stability
        elif hull_distance < 200:
            color_code = 0.71

        # Poor chance of stability
        else:
            color_code = 0.92

        ax.add_patch(plt.Rectangle((i + 0.1, 0), height=hull_distance,
                                   width=0.8, linewidth=0,
                                   facecolor=plt.cm.jet(color_code)))
        i += 1

    ax.set_xticks([x + 0.5 for x in range(len(hull_distances))])
    ax.set_xticklabels(x_ticklabels, family='serif', size=20, rotation=60)
    ax.set_yticklabels(ax.get_yticks(), family='serif', size=20)
    ax.set_ylabel(r'$\mathrm{E_F\/(meV/atom)}$', size=40)

    plt.savefig('stability_plot.{}'.format(fmt), transparent=True)


def all_converged(two_d, three_d):
    """
    True if all two_d directories and three_d competing_phases have
    converged, otherwise false.
    """

    finished_2d, finished_3d = [], []
    for directory in two_d:
        if is_converged(directory):
            finished_2d.append(directory)
    for directory in three_d:
        if is_converged('all_competitors/{}'.format(directory[0])):
            finished_3d.append(directory[0])
    if len(finished_2d + finished_3d) == len(
            two_d + three_d):
        all_converged = True
    else:
        all_converged = False
    return all_converged
