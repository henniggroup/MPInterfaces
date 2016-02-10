import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from pymatgen.io.vasp.outputs import Vasprun


def plot_friction_surface(directory):
    """
    Collect the energies along a 10x10 grid of static energy
    calculations to plot the Gamma surface between two layers of the 2D
    material.
    """

    ax = plt.figure(figsize=(10, 10)).gca()

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)

    ENERGY_ARRAY = []

    X_VALUES = range(10)
    Y_VALUES = range(10)

    for x in X_VALUES:
        ENERGY_ARRAY.append([])
        for y in Y_VALUES:
            dir = '{}x{}'.format(x, y)
            os.chdir(dir)
            ENERGY_ARRAY[x].append(Vasprun('vasprun.xml').final_energy)
            os.chdir('../')

    minima = []
    maxima = []
    for x in X_VALUES:
        minima.append(min(ENERGY_ARRAY[x]))
        maxima.append(max(ENERGY_ARRAY[x]))

    abs_minimum = min(minima)
    abs_maximum = max(maxima)

    for x in X_VALUES:
        for y in Y_VALUES:

            # Plot all energies relative to the global minimum.
            scaled_energy = ENERGY_ARRAY[x][y] - abs_minimum
            color_code = scaled_energy / (abs_maximum - abs_minimum)

            ax.add_patch(plt.Rectangle((x, y), width=1, height=1,
                                       facecolor=plt.cm.coolwarm(color_code),
                                       linewidth=0))
    # Get rid of annoying ticks.
    ax.axes.get_yaxis().set_ticks([])
    ax.axes.get_xaxis().set_ticks([])

    plt.savefig('{}.pdf'.format(os.getcwd().split('/')[-1]),
                transparent='True')
