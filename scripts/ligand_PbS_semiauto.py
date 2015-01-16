"""
This script creates the four input files for an interface relaxation job
pending setting of exact sd flags
"""
import shutil as shu

from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from pymatgen.core import Structure, Molecule


from mpinterfaces.interface import Interface, Ligand

        


#create ligand, interface and slab from the starting POSCARs


strt= Structure.from_file("POSCAR_PbS_bulk_with_vdw")  #using POSCAR of vdW relaxed PbS
mol_struct= Structure.from_file("POSCAR_DMF")    #using POSCAR of vdW relaxed PbS
mol= Molecule(mol_struct.species, mol_struct.cart_coords)
DMF= Ligand([mol])   #create Ligand DMF 
supercell = [1,1,1]
# slab thickness and vacuum set manual for now to converged values, surface coverage fixed at 0.014 ligand/sq.Angstrom 
#for consistency, best ligand spacing at the coverage 
min_thick= 19
min_vac= 12
surface_coverage= 0.014
#hkl of facet to reproduce
hkl= [1,0,0]
# specify the species on slab to adsorb over 
slab_species= 'Pb'
# specify the species onb ligand serving as the bridge atom 
adatom_on_ligand= 'O' 
#initial adsorption distance in angstrom
ads_distance = 3.0 
# create the interface
iface = Interface(strt, hkl=hkl, min_thick=min_thick, min_vac=min_vac,
                      supercell=supercell, surface_coverage=surface_coverage,
                      ligand=DMF, displacement=ads_distance, adatom_on_lig= adatom_on_ligand, adsorb_on_species= slab_species, primitive= True)
iface.create_interface()
iface.sort()
#separate slab 
iface_slab = iface.slab
iface_slab.sort()
#set selective dynamics flags as required
true_site= [1, 1, 1]
false_site= [0, 0, 0]
sd_flag_iface= []
j= 0
sd_flag_slab= []
for i in iface.sites:
	sd_flag_iface.append(false_site)
for i in iface_slab.sites:
        sd_flag_slab.append(false_site)

#POSCAR and POTCAR construction, pending setting of exact flags for POSCAR
interface_poscar = Poscar(iface, selective_dynamics= sd_flag_iface)
interface_potcar= Potcar(interface_poscar.site_symbols)
slab_poscar = Poscar(iface_slab, selective_dynamics= sd_flag_slab)
slab_potcar= Potcar(slab_poscar.site_symbols)

#write the files in appropriate directories 
interface_poscar.write_file("./Interface/POSCAR_PbS100DMF_interface_with_sd.vasp")
slab_poscar.write_file("./Slab/POSCAR_PbS100_slab_with_sd.vasp")
interface_potcar.write_file("./Interface_with_vdw/POTCAR")
slab_potcar.write_file("./Slab_with_vdw/POTCAR")

#set the common INCAR file and KPOINTS
incar_dict = {
                 'SYSTEM': 'ligand_PbS', 
                 'ENCUT': 600, 
                 'ISIF': 2, 
                 'IBRION': 2, 
                 'ALGO': 'Normal', 
                 'ISMEAR': 1, 
                 'ISPIN': 1, 
                 'EDIFF': 1e-06, 
                 'EDIFFG': -0.005, 
                 'NPAR': 8, 
                 'SIGMA': 0.1, 
                 'PREC': 'Accurate',
		 'IVDW': 2,
		 'NSW': 1000
    }

incar = Incar.from_dict(incar_dict)
kpoints = Kpoints.monkhorst_automatic(kpts= (8, 8, 1), shift= (0, 0, 0))

#write the files in appropriate directory
incar.write_file("./Interface_with_vdw/INCAR")
incar.write_file("./Slab_with_vdw/INCAR")
kpoints.write_file("./Interface_with_vdw/KPOINTS")
kpoints.write_file("./Slab_with_vdw/KPOINTS")
shu.copy("./submit_job", "./Interface_with_vdw/")
shu.copy("./submit_job", "./Slab_with_vdw")
