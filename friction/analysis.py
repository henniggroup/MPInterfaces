import os

import matplotlib.pyplot as plt

from pymatgen.io.vasp.outputs import Vasprun


def plot_friction_surface(directory):

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
    print abs_minimum, abs_maximum
    print ENERGY_ARRAY[9][9]

    for x in X_VALUES:
        for y in Y_VALUES:
            scaled_energy = ENERGY_ARRAY[x][y] - abs_minimum
            color_code = scaled_energy / (abs_maximum - abs_minimum)

            ax.add_patch(plt.Rectangle((x, y), width=1, height=1,
                                       facecolor=plt.cm.coolwarm(color_code),
                                       linewidth=0))

    ax.axes.get_yaxis().set_ticks([])
    ax.axes.get_xaxis().set_ticks([])

    plt.savefig('{}.pdf'.format(os.getcwd().split('/')[-1]),
                transparent='True')
