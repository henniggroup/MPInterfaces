import os

from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.matproj.rest import MPRester

from monty.serialization import loadfn

import twod_materials.standard as st

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


MPR = MPRester(
    loadfn(os.path.join(os.path.expanduser('~'), 'dbauth.yaml'))['mp_api']
    )


def get_competing_species(directories):

    total_competing_species = []

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
        competing_species = [
            (entry.composition.reduced_formula,
             entry.entry_id) for entry in decomp[0]
            ]

        # Keep a running list of all unique competing species, since in
        # high throughput 2D searches there is usually some overlap in
        # competing species for different materials.
        for specie in competing_species:
            if specie not in total_competing_species:
                total_competing_species.append(specie)
        os.chdir('../')

    return total_competing_species


def get_hull_distances(directories):

    hull_distances = {}
    finished_competitors = {}

    # Determine which competing species have been relaxed in the current
    # framework and store them in a dictionary ({formula: entry}).
    if os.path.isdir('all_competitors'):
        os.chdir('all_competitors')
        for comp_dir in [
            dir for dir in os.listdir(os.getcwd()) if os.path.isdir(dir) and
            st.is_converged(dir)
                ]:
            os.chdir(comp_dir)
            composition = Structure.from_file('POSCAR').composition
            energy = Vasprun('vasprun.xml').final_energy
            finished_competitors[comp_dir] = ComputedEntry(composition, energy)
            os.chdir('../')
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

        # If the energies of competing species have been calculated in
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

    return hull_distances


def plot_hull_distances(hull_distances):

    ax = plt.figure().gca()
    ax.set_ylim(0, 1000)
    ax.set_xlim(0, len(hull_distances))

    x_ticklabels = []
    i = 0
    for compound in hull_distances:
        x_ticklabels.append(compound)
        if hull_distances[compound] < 100:
            color_code = 0.5
        elif hull_distances[compound] < 200:
            color_code = 0.71
        else:
            color_code = 0.92

        ax.add_patch(plt.Rectangle(i, 0), height=hull_distances[compound],
                     width=1, linewidth=0, facecolor=plt.cm.jet(color_code))

    ax.set_xticklabels(x_ticklabels)

    plt.savefig('stability_plot.pdf', transparent=True)
