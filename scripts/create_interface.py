"""
This script demonstrates the usage of the module mpinterfaces/interface.py

Also, demonstrates how to fetch data from the materialsproject database using their API
Note: Before using the script, make sure that you do have a valid api key obtained from the materialsproject website. Use that to set the MAPI_KEY variable below
"""

import sys
from pymatgen.matproj.rest import MPRester
from pymatgen.core import Molecule
from mpinterfaces.interface import Interface, Ligand

MAPI_KEY="dwvz2XCFUEI9fJiR"

def get_struct_from_mp(formula):
    """
    fetches the structure corresponding to the given formula
    from the materialsproject database
    Note: get the api key from materialsproject website
    provide the api key here os set the environment variable "MAPI_KEY"
    Note: for the given formula there are many structures available, this
    function returns the first one of those structures
    """
    with MPRester(MAPI_KEY) as m:
        #criteria = {'elements': {'$in': ['Li', 'Na', 'K'], '$all': ['O']}}
        #props = ['pretty_formula', 'energy']
        #data = m.query(criteria=criteria, properties=props)
        #data = m.get_exp_thermo_data("Fe2O3")
        #structs = m.get_structures("Mn3O4")
        #ntries = m.get_entries("TiO2")
        #entry = m.get_exp_entry("Fe2O3")
        data = m.get_data(formula)
        print "\nnumber of structures matching the chemical formula "+formula+" = ", len(data)
        #srtrs = {}
        #get info for each structure
        for d in data:
            x = {}
            x['material_id'] = str(d['material_id'])
            structure = m.get_structure_by_material_id(x['material_id'])
            return structure
            #
            #x['spacegroup'] = str(d['spacegroup']['symbol'])
            #x['formation_energy_per_atom'] = d['formation_energy_per_atom']
            #x['band_gap'] = d['band_gap']
            #srtrs[str(d['full_formula'])] = x
            #print "Unit cell vol = {}".format(structure.volume), d['volume']
            #print the info
            #for k,v in srtrs.iteritems():
            #    print "\n", k, " : \n"
            #    for k1,v1 in srtrs[k].iteritems():
            #        print k1, " : ", v1
            #Dos for material id
            #  dos = m.get_dos_by_material_id("mp-1234")
            #Bandstructure for material id
            #  bandstructure = m.get_bandstructure_by_material_id("mp-1234")


                                                                                       
if __name__=='__main__':
    # PbS 100 surface with single hydrazine as ligand
    strt= Structure.from_file("POSCAR_PbS")  #provide POSCAR of building block POSCAR_bulk and POSCAR_molecule here
    mol_struct= Structure.from_file("POSCAR_Hydrazine")
    mol= Molecule(mol_struct.species, mol_struct.cart_coords)
    hydrazine= Ligand([mol])

    #intital supercell, this wont be the final supercell if surface coverage is specified
    supercell = [1,1,1]

    #miller index
    hkl = [1,0,0]
    
    #minimum slab thickness in Angstroms
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
    iface = Interface(strt, hkl=hkl, min_thick=min_thick, min_vac=min_vac,
                      supercell=supercell, surface_coverage=0.01,
                      ligand=hydrazine, displacement=displacement, adatom_on_lig='N',
                      adsorb_on_species= 'Pb', primitive= False)
    iface.create_interface()
    iface.sort()
    #extract bare slab
    iface_slab = iface.slab
    iface_slab.sort()
    #set selective dynamics flags as required
    true_site= [1, 1, 1]
    false_site= [0, 0, 0]
    sd_flag_iface= []
    sd_flag_slab= []
    #selective dynamics flags for the interface
    for i in iface.sites:
    	sd_flag_iface.append(false_site)
    #selective dynamics flags for the bare slab        
    for i in iface_slab.sites:
        sd_flag_slab.append(false_site)
    interface_poscar = Poscar(iface, selective_dynamics= sd_flag_iface)
    slab_poscar = Poscar(iface_slab, selective_dynamics= sd_flag_slab)
    #poscars without selective dynamics flag
    iface.to('poscar', 'POSCAR_interface.vasp')
    iface_slab.to('poscar', 'POSCAR_slab.vasp')
    #poscars with selective dynamics flag    
    interface_poscar.write_file("POSCAR_interface_with_sd.vasp")
    slab_poscar.write_file("POSCAR_slab_with_sd.vasp") 
