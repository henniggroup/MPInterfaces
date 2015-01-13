"""
Goes through a typical volume data stored in Directory VOLUME and obtains the EOS plot of the E(V) data according to Birch 
Murnaghan

TODO: Incorporate borg assimilate for the large directory tree structure
	formatting with matplotlib optional eos.plot gives a .png image 

"""

import numpy as np
import os
from pymatgen.io.vaspio import Vasprun
from pymatgen.core.structure import Structure, Lattice
#import matplotlib.pyplot as plot
from ase.utils.eos import EquationOfState

#data_file = file("DATA.dat", 'w') #optional data file to write to 
energy_data= []
volume_data= []
for i in ["0_97", "0_98", "0_99", "1_0", "1_01", "1_02", "1_03"]: #replace with sequence of volume directory names
	name = str(i)
	path_name = "./VOLUME/" + name
	os.chdir(path_name)
	obj = Vasprun(filename= "vasprun.xml")
	os.chdir("../../")
	struct= obj.structures
	energy= obj.final_energy
	volume= struct[0].lattice.volume
	energy_data.append(energy)
	volume_data.append(volume) 
#volumes= np.array(volume_data)
#energies= np.array(energy_data)
print "plotting: Energy vs Volume", energy_data, volume_data
eos = EquationOfState(volumes= volume_data, energies= energy_data, eos='murnaghan')
v0, e0, B = eos.fit()
#the ASE units for the bulk modulus is eV/Angstrom^3 
print 'optimum volume, energy(eV) and bulk modulus(eV/Angstrom^3)', v0, e0, B
    #plot
eos.plot(filename= "eos_fit_plot")    
    	#plot.plot(volumes,energies,'ro')
    	#x = np.linspace(min(eos.v), max(eos.v), 100)
    	#y = eval(eos.eos_string)(x, eos.eos_parameters[0], eos.eos_parameters[1],
        #                     eos.eos_parameters[2], eos.eos_parameters[3] )
    	#plot.plot(x, y, label='fit')
    	#plot.xlabel('Volume ($\AA^3$)')
    	#plot.ylabel('Energy (eV)')
    	#plot.legend(loc='best')
    	#plot.savefig('eos.png')
    #show()

