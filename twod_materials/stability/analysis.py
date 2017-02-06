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

import twod_materials


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')
PACKAGE_PATH = '/'.join(PACKAGE_PATH.split('/')[:-2])

try:
    config_vars = loadfn(os.path.join(os.path.expanduser('~'), 'config.yaml'))
except:
    print('WARNING: No config.yaml file was found. please configure the '\
    'config.yaml and put it in your home directory.')
    # Still set them for testing purposes.
    config_vars = loadfn(os.path.join(PACKAGE_PATH, 'config.yaml'))
if 'MP_API' in os.environ:
    MPR = MPRester(os.environ['MP_API'])
else:
    MPR = MPRester(config_vars['mp_api'])


def get_competing_phases():
    """
    Collect the species to which the material might decompose to.

    Returns:
        A list of phases as tuples formatted as
        [(formula_1, Materials_Project_ID_1),
        (formula_2, Materials_Project_ID_2), ...]
    """

    total_competing_phases = []

    composition = Structure.from_file('POSCAR').composition
    try:
        energy = Vasprun('vasprun.xml').final_energy
    except:
        energy = 100  # The function can work without a vasprun.xml
    entries = MPR.get_entries_in_chemsys(
        [elt.symbol for elt in composition]
        )
    my_entry = ComputedEntry(composition, energy)
    entries.append(my_entry)

    pda = PDAnalyzer(PhaseDiagram(entries))
    decomp = pda.get_decomp_and_e_above_hull(my_entry, allow_negative=True)
    competing_phases = [
        (entry.composition.reduced_formula,
         entry.entry_id) for entry in decomp[0]
        ]

    return competing_phases


def get_hull_distance(competing_phase_directory='../competing_phases'):
    """
    Calculate the material's distance to the thermodynamic hull,
    based on species in the Materials Project database.

    Args:
        competing_phase_directory (str): absolute or relative path
            to the location where your competing phases have been
            relaxed. The default expectation is that they are stored
            in a directory named 'competing_phases' at the same level
            as your material's relaxation directory.
    Returns:
        float. distance (eV/atom) between the material and the
            hull.
    """

    finished_competitors = {}
    original_directory = os.getcwd()
    # Determine which competing phases have been relaxed in the current
    # framework and store them in a dictionary ({formula: entry}).
    if os.path.isdir(competing_phase_directory):
        os.chdir(competing_phase_directory)
        for comp_dir in [
            dir for dir in os.listdir(os.getcwd()) if os.path.isdir(dir) and
            is_converged(dir)
                ]:
            vasprun = Vasprun('{}/vasprun.xml'.format(comp_dir))
            composition = vasprun.final_structure.composition
            energy = vasprun.final_energy
            finished_competitors[comp_dir] = ComputedEntry(composition, energy)
        os.chdir(original_directory)
    else:
        raise ValueError('Competing phase directory does not exist.')

    composition = Structure.from_file('POSCAR').composition
    try:
        energy = Vasprun('vasprun.xml').final_energy
    except:
        raise ValueError('This directory does not have a converged vasprun.xml')
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

    return decomp[1]


def plot_hull_distances(hull_distances, fmt='pdf'):
    """
    Create a bar graph of the formation energies of several 2D materials.

    Args:
        hull_distances (dict): follow the format:
            {reduced_formula: hull_distance (in eV/atom)}
        fmt (str): matplotlib format style. Check the matplotlib
            docs for options.
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
