"""
Relaxes 2D materials in all subdirectories of the current working
directory, along with their most stable competing species. At a
specified INTERVAL, checks if all relaxations have converged. Once all
are converged, calculates and plots the formation energies of all 2D
materials as stability_plot.pdf.
"""
from __future__ import print_function, division, unicode_literals

import os

import time

from twod_materials.utils import is_converged
from twod_materials.electronic_structure.startup import (
    run_linemode_calculation
    )
from twod_materials.electronic_structure.analysis import (
    plot_normal_band_structure
    )

INTERVAL = 360  # Seconds between convergence checks

directories = [dir for dir in os.listdir(os.getcwd()) if os.path.isdir(dir)
               and dir not in ['all_competitors']]

if __name__ == '__main__':

    for directory in directories:
        os.chdir(directory)
        run_linemode_calculation()
        os.chdir('../')

    loop = True
    while loop:
        print('>> Checking convergence')
        finished = []

        for directory in directories:
            if is_converged('{}/pbe_bands'.format(directory)):
                finished.append(directory)

        if len(finished) == len(directories):
            print('>> Plotting band structures')
            for directory in finished:
                os.chdir('{}/pbe_bands'.format(directory))
                plot_normal_band_structure()
                os.chdir('../../')
            loop = False
        else:
            print('>> Not all directories converged ({}/{})'.format(
                len(finished), len(directories)))

            time.sleep(INTERVAL)
