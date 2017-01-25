from __future__ import print_function, division, unicode_literals

from twod_materials.friction.startup import run_friction_calculations
from twod_materials.friction.analysis import plot_gamma_surface

import os

import time

from twod_materials.utils import is_converged


INTERVAL = 360  # Seconds between convergence checks

directories = [dir for dir in os.listdir(os.getcwd()) if os.path.isdir(dir)
               and dir not in ['all_competitors']]

if __name__ == '__main__':

    converged = {}

    for directory in directories:
        os.chdir(directory)
        converged[directory] = False
        run_friction_calculations()
        os.chdir('../')

    loop = True

    while loop:
        time.sleep(INTERVAL)

        loop = False

        for directory in directories:
            os.chdir(directory)
            converged[directory] = True
            for subdirectory in [
                    dir for dir in os.listdir(os.getcwd())
                    if os.path.isdir(dir)]:
                if not is_converged(subdirectory):
                    converged[directory] = False
                    break

            if not converged[directory]:
                print('>> Not all directories converged')
                loop = True

            print('>> Plotting gamma surface for {}'.format(directory))
            plot_gamma_surface()

            os.chdir('../')
