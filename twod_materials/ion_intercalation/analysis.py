import os

from twod_materials.utils import is_converged

from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.io.vasp.outputs import Vasprun

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import operator


def plot_lithium_hull():
    """
    Plots the phase diagram between the pure material and pure Li,
    Connecting the points on the convex hull of the phase diagram.
    """

    # Calculated with the relax_3d() function
    # in twod_materials.stability.startup
    Li_ev_fu = -1.7540797

    energy = Vasprun('vasprun.xml').final_energy
    composition = Structure.from_file('POSCAR').composition
    twod_material = composition.reduced_formula
    twod_eV_fu = energy / composition.get_reduced_composition_and_factor()[1]

    if os.path.isdir('lithium'):
        os.chdir('lithium')

        hull = {0: 0, 1: 0}  # Set both end points to zero E_formation

        for directory in [
                dir for dir in os.listdir(os.getcwd()) if os.path.isdir(dir)]:
            if is_converged(directory):
                os.chdir(directory)
                energy = Vasprun('vasprun.xml').final_energy
                composition = Structure.from_file('POSCAR').composition
                Li_fraction = composition.get_atomic_fraction('Li')

                no_Li_comp_dict = composition.as_dict()
                no_Li_comp_dict.update({'Li': 0})
                no_Li_comp = Composition.from_dict(no_Li_comp_dict)

                n_twod_fu = no_Li_comp.get_reduced_composition_and_factor()[1]

                hull[Li_fraction] = (
                    (energy - composition['Li'] * Li_ev_fu
                     - twod_eV_fu * n_twod_fu)
                    / composition.num_atoms
                )
                os.chdir('../')
        os.chdir('../')

    sorted_hull = sorted(hull.items(), key=operator.itemgetter(0))

    # Determine which compositions are on the convex hull.
    convex_points = [(0, 0)]
    concave_points = []
    n = 1
    for i in range(1, len(sorted_hull) - 1):
        in_slope = (
            (sorted_hull[i][1] - convex_points[n - 1][1])
            / (sorted_hull[i][0] - convex_points[n - 1][0])
        )

        out_slope = (
            (sorted_hull[i][1] - sorted_hull[i + 1][1])
            / (sorted_hull[i][0] - sorted_hull[i + 1][0])
        )

        if out_slope >= in_slope:
            convex_points.append(sorted_hull[i])
            n += 1
        else:
            concave_points.append(sorted_hull[i])
    convex_points.append((1, 0))

    convex_li_fractions = [tup[0] for tup in convex_points]
    convex_formation_energies = [tup[1] for tup in convex_points]
    concave_li_fractions = [tup[0] for tup in concave_points]
    concave_formation_energies = [tup[1] for tup in concave_points]

    ax = plt.figure(figsize=(10, 10)).gca()

    ax.plot([0, 1], [0, 0], 'k--')
    ax.plot(convex_li_fractions, convex_formation_energies, 'b-', marker='o',
            markersize=12)
    ax.plot(concave_li_fractions, concave_formation_energies, 'r', marker='o',
            linewidth=0, markersize=12)

    ax.text(0, 0.002, twod_material, family='serif', size=20)
    ax.text(0.99, 0.002, 'Li', family='serif', size=20,
            horizontalalignment='right')

    ax.set_xlabel(r'$\mathrm{at\/\%\/Li}$', size=20)
    ax.set_ylabel(r'$\mathrm{E_F\/(eV/atom)}$', size=20)

    plt.savefig('lithium_hull.eps')
