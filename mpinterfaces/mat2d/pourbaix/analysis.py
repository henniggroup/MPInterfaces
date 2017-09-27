"""
Create pourbaix diagrams for 2D materials against their ions in solution using
the scheme outlined in PHYSICAL REVIEW B 85, 235438 (2012)
"""
from __future__ import print_function, division, unicode_literals

import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from monty.serialization import loadfn

from pymatgen import Element
from pymatgen.analysis.pourbaix.analyzer import PourbaixAnalyzer
from pymatgen.analysis.pourbaix.entry import PourbaixEntry, IonEntry
from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
from pymatgen.analysis.pourbaix.plotter import PourbaixPlotter
from pymatgen.core.ion import Ion
from pymatgen.core.structure import Structure
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.io.vasp.outputs import Vasprun

import seaborn as sb


__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


PACKAGE_PATH = os.path.dirname(__file__)
ION_CORRECTIONS = loadfn(os.path.join(PACKAGE_PATH, 'ion_corrections.yaml'))
ION_FORMATION_ENERGIES = loadfn(os.path.join(PACKAGE_PATH,
                                             'ion_formation_energies.yaml'))
CHEMICAL_POTENTIALS = loadfn(os.path.join(PACKAGE_PATH,
                                          'chemical_potentials.yaml'))


def contains_entry(entry_list, entry):
    """
    Function for filtering duplicate entries in a list.

    Args:
        entry_list (list): List of Pymatgen ComputedEntry objects.
        entry (ComputedEntry): the Pymatgen ComputedEntry object of
            interest.

    Returns:
        bool
    """
    for ent in entry_list:
        if (ent.entry_id == entry.entry_id or
                (abs(entry.energy_per_atom - ent.energy_per_atom) < 1e-6 and
                     (entry.composition.reduced_formula ==
                      ent.composition.reduced_formula))):
            return True


def plot_pourbaix_diagram(metastability=0.0, ion_concentration=1e-6, fmt='pdf'):
    """
    Creates a Pourbaix diagram for the material in the cwd.

    Args:
        metastability (float): desired metastable tolerance energy
            (meV/atom). <~50 is generally a sensible range to use.
        ion_concentration (float): in mol/kg. Sensible values are
            generally between 1e-8 and 1.
        fmt (str): matplotlib format style. Check the matplotlib
            docs for options.
    """

    # Create a ComputedEntry object for the 2D material.
    composition = Structure.from_file('POSCAR').composition
    energy = Vasprun('vasprun.xml').final_energy

    cmpd = ComputedEntry(composition, energy)

    # Define the chemsys that describes the 2D compound.
    chemsys = ['O', 'H'] + [elt.symbol for elt in composition.elements
                            if elt.symbol not in ['O', 'H']]

    # Pick out the ions pertaining to the 2D compound.
    ion_dict = dict()
    for elt in chemsys:
        if elt not in ['O', 'H'] and ION_FORMATION_ENERGIES[elt]:
            ion_dict.update(ION_FORMATION_ENERGIES[elt])

    elements = [Element(elt) for elt in chemsys if elt not in ['O', 'H']]

    # Add "correction" for metastability
    cmpd.correction -= float(cmpd.composition.num_atoms)\
        * float(metastability) / 1000.0

    # Calculate formation energy of the compound from its end
    # members
    form_energy = cmpd.energy
    for elt in composition.as_dict():
        form_energy -= CHEMICAL_POTENTIALS[elt] * cmpd.composition[elt]

    # Convert the compound entry to a pourbaix entry.
    # Default concentration for solid entries = 1
    pbx_cmpd = PourbaixEntry(cmpd)
    pbx_cmpd.g0_replace(form_energy)
    pbx_cmpd.reduced_entry()

    # Add corrected ionic entries to the pourbaix diagram
    # dft corrections for experimental ionic energies:
    # Persson et.al PHYSICAL REVIEW B 85, 235438 (2012)
    pbx_ion_entries = list()

    # Get PourbaixEntry corresponding to each ion.
    # Default concentration for ionic entries = 1e-6
    # ion_energy = ion_exp_energy + ion_correction * factor
    # where factor = fraction of element el in the ionic entry
    # compared to the reference entry
    for elt in elements:
        for key in ion_dict:
            comp = Ion.from_formula(key)
            if comp.composition[elt] != 0:
                factor = comp.composition[elt]
                energy = ion_dict[key]
                pbx_entry_ion = PourbaixEntry(IonEntry(comp, energy))
                pbx_entry_ion.correction = (
                    ION_CORRECTIONS[elt.symbol] * factor
                )
                pbx_entry_ion.conc = ion_concentration
                pbx_entry_ion.name = key
                pbx_ion_entries.append(pbx_entry_ion)

    # Generate and plot Pourbaix diagram
    # Each bulk solid/ion has a free energy g of the form:
    # g = g0_ref + 0.0591 * log10(conc) - nO * mu_H2O +
    # (nH - 2nO) * pH + phi * (-nH + 2nO + q)

    all_entries = [pbx_cmpd] + pbx_ion_entries

    total = sum([composition[el] for el in elements])
    comp_dict = {el.symbol: composition[el]/total for el in elements}
    pourbaix_diagram = PourbaixDiagram(all_entries, comp_dict=comp_dict)

    plotter = PourbaixPlotter(pourbaix_diagram)

    # Plotting details...
    font = "serif"
    fig = plt.figure(figsize=(14, 9))
    ax1 = fig.gca()
    ax1.set_xlim([0, 14])
    ax1.set_xticklabels([int(t) for t in ax1.get_xticks()], fontname=font,
                        fontsize=18)
    ax1.set_ylim(-2, 2)
    ax1.set_yticklabels(ax1.get_yticks(), fontname=font, fontsize=18)
    ax1.set_xlabel("pH", fontname=font, fontsize=18)
    ax1.set_ylabel("Potential vs. SHE (V)", fontname=font, fontsize=18)

    # Outline water's stability range.
    ax1.plot([0, 14], [0, -0.829], color="gray", linestyle="--", alpha=0.7,
             linewidth=2)
    ax1.plot([0, 14], [1.229, 0.401], color="gray", linestyle="--", alpha=0.7,
             linewidth=2)


    stable_entries = plotter.pourbaix_plot_data(
        limits=[[0, 14], [-2, 2]])[0]

    # Add coloring.
    colors = sb.color_palette("Set2", len(stable_entries))

    i = 0
    for entry in stable_entries:
        col = colors[i]
        i += 1
        vertices = plotter.domain_vertices(entry)
        center_x = sum([v[0] for v in vertices])/len(vertices)
        center_y = sum([v[1] for v in vertices])/len(vertices)
        patch = Polygon(vertices, closed=True, fill=True, facecolor=col,
                        linewidth=2, edgecolor="w")
        ax1.text(center_x, center_y, plotter.print_name(entry),
                 verticalalignment="center", horizontalalignment="center",
                 fontname=font, fontsize=18)
        ax1.add_patch(patch)

    plt.savefig("pourbaix.{}".format(fmt))

    plt.close()
