from __future__ import division, unicode_literals, print_function

"""
Fit energy and volume to an equation of state to get the minimum
volume and the bulk modulus

uses the ase package to do the fitting
"""

import numpy as np
import matplotlib.pyplot as plot

from ase.utils.eos import EquationOfState
from ase.utils.eosase2 import *

from pymatgen.apps.borg.queen import BorgQueen
from mpinterfaces import MPINTVaspDrone

def get_e_v(fname):
    data = np.loadtxt(fname, usecols = (1,3,5,6,7))    
    volumes = data[:,1]
    energies = data[:,4]
    return volumes, energies

def custom_plot(volumes, energies, eos):
    plot.plot(volumes,energies,'ro')
    x = np.linspace(min(eos.v), max(eos.v), 100)
    y = eval(eos.eos_string)(x, eos.eos_parameters[0], eos.eos_parameters[1],
                             eos.eos_parameters[2], eos.eos_parameters[3] )
    plot.plot(x, y, label='fit')
    plot.xlabel('Volume ($\AA^3$)')
    plot.ylabel('Energy (eV)')
    plot.legend(loc='best')
    plot.savefig('eos.png')
    plot.show()

    
if __name__=='__main__':
    volumes, energies = get_e_v('test_bulk.txt')
    #eos = 'sjeos', 'murnaghan', 'birch', 'taylor', 'vinet' etc.
    eos = EquationOfState(volumes, energies, eos='murnaghan')
    v0, e0, B = eos.fit()
    #the ASE units for the bulk modulus is eV/Angstrom^3 
    print('optimum volume, energy and bulk moduls', v0, e0, B)
    #plot
    #eos.plot(filename= "eos_fit")
    custom_plot(volumes, energies, eos)

