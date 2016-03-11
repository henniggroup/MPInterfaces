import os

import numpy as np

import math

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.structure import Structure

import warnings


def plot_gamma_surface():
    """
    Collect the energies along a 10x10 grid of static energy
    calculations to plot the Gamma surface between two layers of the 2D
    material.
    """

    lattice = Structure.from_file('CONTCAR').lattice
    area = np.cross(lattice._matrix[0], lattice._matrix[1])[2]

    ax = plt.figure(figsize=(10, 10)).gca()

    n_divs_x = int(math.ceil(lattice.a * 2.5))
    n_divs_y = int(math.ceil(lattice.b * 2.5))

    ax.set_xlim(0, n_divs_x)
    ax.set_ylim(0, n_divs_y)

    ENERGY_ARRAY = []

    X_VALUES = range(n_divs_x)
    Y_VALUES = range(n_divs_y)

    os.chdir('friction')

    not_converged = []
    for x in X_VALUES:
        ENERGY_ARRAY.append([])
        for y in Y_VALUES:
            dir = '{}x{}'.format(x, y)
            os.chdir(dir)
            try:
                energy = Vasprun('vasprun.xml').final_energy / area
                ENERGY_ARRAY[x].append(energy)
            except:
                not_converged.append('{}x{}'.format(x, y))
                ENERGY_ARRAY[x].append(0)
            os.chdir('../')
        ENERGY_ARRAY[x].append(ENERGY_ARRAY[x][0])

    os.chdir('../')

    ENERGY_ARRAY.append([])
    ENERGY_ARRAY[n_divs_x] = ENERGY_ARRAY[0]

    if not_converged:
        warnings.warn('{} did not converge.'.format(not_converged))
        for coords in not_converged:
            ENERGY_ARRAY[
                int(coords.split('x')[0])][int(coords.split('x')[1])] = energy

    minima = []
    maxima = []
    for x in X_VALUES:
        minima.append(min(ENERGY_ARRAY[x]))
        maxima.append(max(ENERGY_ARRAY[x]))

    abs_minimum = min(minima)
    abs_maximum = max(maxima)

    for x in range(n_divs_x + 1):
        for y in range(n_divs_y + 1):

            # Plot all energies relative to the global minimum.
            scaled_energy = ENERGY_ARRAY[x][y] - abs_minimum
            if '{}x{}'.format(x, y) in not_converged:
                color_code = 'w'
            else:
                color_code = plt.cm.jet(
                    scaled_energy / (abs_maximum - abs_minimum)
                    )

            # Label the energy if it's a saddle point or local maximum.
            try:
                if (
                    (ENERGY_ARRAY[x-1][y] > ENERGY_ARRAY[x][y] <
                     ENERGY_ARRAY[x+1][y] and ENERGY_ARRAY[x][y-1] <
                     ENERGY_ARRAY[x][y] > ENERGY_ARRAY[x][y+1])
                    or
                    (ENERGY_ARRAY[x-1][y] < ENERGY_ARRAY[x][y] >
                     ENERGY_ARRAY[x+1][y] and ENERGY_ARRAY[x][y-1] >
                     ENERGY_ARRAY[x][y] < ENERGY_ARRAY[x][y+1])
                    or
                    (ENERGY_ARRAY[x-1][y-1] > ENERGY_ARRAY[x][y] <
                     ENERGY_ARRAY[x+1][y+1] and ENERGY_ARRAY[x+1][y-1] <
                     ENERGY_ARRAY[x][y] > ENERGY_ARRAY[x-1][y+1])
                    or
                    (ENERGY_ARRAY[x-1][y] < ENERGY_ARRAY[x][y] >
                     ENERGY_ARRAY[x+1][y] and ENERGY_ARRAY[x][y-1] >
                     ENERGY_ARRAY[x][y] < ENERGY_ARRAY[x][y+1])
                    or
                    (ENERGY_ARRAY[x][y] > max([ENERGY_ARRAY[x-1][y-1],
                                               ENERGY_ARRAY[x+1][y+1],
                                               ENERGY_ARRAY[x-1][y+1],
                                               ENERGY_ARRAY[x+1][y-1],
                                               ENERGY_ARRAY[x-1][y],
                                               ENERGY_ARRAY[x+1][y],
                                               ENERGY_ARRAY[x][y-1],
                                               ENERGY_ARRAY[x][y+1]]))
                        ):
                    ax.text(x+0.45, y+0.45, round(scaled_energy*1000, 2),
                            family='serif', size=20,
                            horizontalalignment='center',
                            verticalalignment='center')
            except IndexError:
                pass

            ax.add_patch(plt.Rectangle((x, y), width=1, height=1,
                                       facecolor=color_code,
                                       linewidth=0))
    # Get rid of annoying ticks.
    ax.axes.get_yaxis().set_ticks([])
    ax.axes.get_xaxis().set_ticks([])

    plt.savefig('gamma_surface.pdf', transparent=True)
