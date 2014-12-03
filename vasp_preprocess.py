import os
import sys
import shutil as shu
from pymatgen.io.vaspio.vasp_input import Incar, Potcar, Kpoints
from pymatgen.core.structure import Structure
import subprocess as sp

#dictionary of defaults for molecule ENCUT convergence
default_incar = {'SYMPREC': 1e-06, 
                 'SYSTEM': 'molecule', 
                 'ENCUT': 500, 
                 'ISIF': 2, 
                 'IBRION': 2, 
                 'ALGO': 'Normal', 
                 'ISMEAR': 1, 
                 'ISPIN': 1, 
                 'EDIFF': 1e-06, 
                 'EDIFFG': -0.001, 
                 'NPAR': 8, 
                 'SIGMA': 0.1, 
                 'PREC': 'Accurate'}

incar = Incar.from_dict(default_incar)
poscar = Structure.from_file("POSCAR", sort=True)
symbols = poscar.symbol_set
potcar = Potcar(symbols)
kpoints = Kpoints()

#create INCAR loop according to ENCUT
#assign ENCUT value
for E in range(400,850,50):
	sp.check_call("rm -rf "+str(E), shell = True)	
	os.mkdir(str(E))
	os.chdir(str(E))
	incar['ENCUT'] = E
	incar.write_file("INCAR")
	poscar.write_file("POSCAR")
	potcar.write_file("POTCAR")
	kpoints.write_file("KPOINTS")
        shu.copy("../submit_job",".")
	ret_val = sp.check_call("qsub submit_job", shell=True)
        if ret_val != 0:
                sys.exit()
	os.chdir('..')




