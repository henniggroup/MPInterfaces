"""
Fit energy and volume to an equation of state to get the minimum
volume and the bulk modulus

uses the ase package to do the fitting
"""
import numpy as np
import matplotlib.pyplot as plot

from ase.utils.eos import EquationOfState
from ase.utils.eosase2 import *

if __name__=='__main__':
    #volumes = np.loadtxt('E-V')[:,0]
    #energies = np.loadtxt('E-V')[:,1]
    volumes = np.array([13.72, 14.83, 16.0, 17.23, 18.52])
    energies = np.array([-56.29, -56.41, -56.46, -56.46, -56.42])
    #eos = 'sjeos', 'murnaghan', 'birch', 'taylor', 'vinet' etc.
    eos = EquationOfState(volumes, energies, eos='murnaghan')
    v0, e0, B = eos.fit()
    #the ASE units for the bulk modulus is eV/Angstrom^3 
    print 'optimum volume, energy and bulk moduls', v0, e0, B
    #plot
    #eos.plot()    
    plot.plot(volumes,energies,'ro')
    x = np.linspace(min(eos.v), max(eos.v), 100)
    y = eval(eos.eos_string)(x, eos.eos_parameters[0], eos.eos_parameters[1],
                             eos.eos_parameters[2], eos.eos_parameters[3] )
    plot.plot(x, y, label='fit')
    plot.xlabel('Volume ($\AA^3$)')
    plot.ylabel('Energy (eV)')
    plot.legend(loc='best')
    plot.savefig('eos.png')
    #show()

