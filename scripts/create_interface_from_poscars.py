"""
This script demonstrates the usage of the module mpinterfaces/interface.py
to make an interface out of an input POSCAR_bulk and POSCAR_molecule

(useful for quick manual tests)
"""

from pymatgen.core import Molecule, Structure
from mpinterfaces.interface import Interface, Ligand
from pymatgen.io.vaspio import Poscar


                                                                                       
if __name__=='__main__':
    # PbS 100 surface with single hydrazine as ligand
    strt= Structure.from_file("POSCAR_PbS")  #provide POSCAR of building block POSCAR_bulk and POSCAR_molecule here
    mol_struct= Structure.from_file("POSCAR_Hydrazine")
    mol= Molecule(mol_struct.species, mol_struct.cart_coords)
    hydrazine= Ligand([mol])
    supercell = [1,1,1]

    #miller index: specify which surface of PbS to create
    hkl = [1,0,0]
    
    #minimum slab thickness in Angstroms (corresponds to 8 layers) 
    min_thick = 19
    
    #minimum vacuum thickness in Angstroms
    #mind: the ligand will be placed in this vacuum, so the
    #final effective vacuum space will be smaller than this
    min_vac = 12
    
    # surface coverage in the units of lig/ang^2
    #mind: exact coverage as provided cannot be guaranteed, the slab will be constructed
    #with a coverage value thats close to the requested one
    #note: maximum supercell size possible is 10 x 10
    #note: 1 lig/nm^2 = 0.01 lig/ang^2    
    surface_coverage = 0.01
    
    #atom on the slab surface on which the ligand will be attached,
    #no need to specify if the slab is made of only a single species
    adsorb_on_species = 'Pb'
    
    #atom on ligand that will be attached to the slab surface
    adatom_on_lig='N'
    
    #ligand displacement from the slab surface along the surface normal
    #i.e adatom_on_lig will be displced by this amount from the adsorb_on_species atom
    #on the slab
    #in Angstrom
    displacement = 3.0

    #
    #here we create the interface
    #
    iface = Interface(strt, hkl=[1,0,0], min_thick=min_thick, min_vac=min_vac,
                      supercell=supercell, surface_coverage=0.01,
                      ligand=hydrazine, displacement=displacement, adatom_on_lig='N', adsorb_on_species= 'Pb', primitive= False)
    iface.create_interface()
    iface.sort()
    #extract bare slab
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
    interface_poscar = Poscar(iface, selective_dynamics= sd_flag_iface)
    slab_poscar = Poscar(iface_slab, selective_dynamics= sd_flag_slab)
    #for poscars that require manual addition of SD_flags
    iface.to('poscar', 'POSCAR_interface.vasp')
    iface_slab.to('poscar', 'POSCAR_slab.vasp')
    #for poscars with sd_flags set
    interface_poscar.write_file("POSCAR_interface_with_sd.vasp")
    slab_poscar.write_file("POSCAR_slab_with_sd.vasp") 
