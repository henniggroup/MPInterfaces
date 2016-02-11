"""
Relaxes 2D materials in all subdirectories of the current working
directory, along with their most stable competing species. At a
specified INTERVAL, checks if all relaxations have converged. Once all
are converged, calculates and plots the formation energies of all 2D
materials as stability_plot.pdf.
"""

import os

import time

from twod_materials.utils import is_converged
from twod_materials.stability.startup import relax, relax_competing_species
from twod_materials.stability.analysis import get_competing_species
from twod_materials.stability.analysis import get_hull_distances
from twod_materials.stability.analysis import plot_hull_distances


INTERVAL = 360  # Seconds between convergence checks

directories = [dir for dir in os.listdir(os.getcwd()) if os.path.isdir(dir)
               and dir not in ['all_competitors']]

if __name__ == '__main__':

    relax(directories)
    competing_species = get_competing_species(directories)
    relax_competing_species(competing_species)

    loop = True
    while loop:
        print '>> Checking convergence'
        finished = []
        for directory in directories:
            if is_converged(directory):
                finished.append(directory)
        for directory in competing_species:
            if is_converged('all_competitors/{}'.format(directory[0])):
                finished.append(directory[0])

        if len(finished) == len(directories + competing_species):
            print '>> Plotting hull_distances'
            plot_hull_distances(get_hull_distances(finished))
            loop = False
        else:
            print '>> Not all directories converged ({}/{})'.format(
                len(finished), len(directories + competing_species))

            time.sleep(INTERVAL)
