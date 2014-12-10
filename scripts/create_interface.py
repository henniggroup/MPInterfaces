"""
This script demonstartes how to fetch data from the materialsproject database using their API
and to create an interface with ligands attached to it
Note: Before using the script, make sure that you do have a valid api key obtained from the materialsproject website
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
    # Cu 111 surface with H2O as ligands

    #h2o ligand
    mol = Molecule(['O','H', 'H'], [ [0,0,0], [0, 0.77, 0.60], [0, -0.77, 0.60]])
    h2o = Ligand([mol])
    
    #Cu bulk structure from materialsproject
    formula = 'Cu' #"Li2O"#"Fe2O3"#"SiO2"#

    #initial structure, must be either a bulk structure or a slab
    strt = get_struct_from_mp(formula)
    
    #intital supercell, this wont be the final supercell if surface coverage is specified
    supercell = [1,1,1]

    #miller index
    hkl = [1,1,1]
    
    #minimum slab thickness in Angstroms
    min_thick = 9
    
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
    #adsorb_on_species = 'Cu'
    
    #atom on ligand that will be attached to the slab surface
    adatom_on_lig='O'
    
    #ligand displacement from the slab surface along the surface normal
    #i.e adatom_on_lig will be displced by this amount from the adsorb_on_species atom
    #on the slab
    #in Angstrom
    displacement = 2.0

    #
    #here we create the interface
    #
    iface = Interface(strt, hkl=[1,1,1], min_thick=min_thick, min_vac=min_vac,
                      supercell=supercell, surface_coverage=0.01,
                      ligand=h2o, displacement=displacement, adatom_on_lig='O')
#    iface = Interface(strt, hkl=hkl, min_thick=min_thick, min_vac=20,
#                      supercell=supercell, surface_coverage=0.01,
#                      ligand=lead_acetate, displacement=displacement, adatom_on_lig='Pb')
    iface.create_interface()
    iface.sort()
    iface.to('poscar', 'POSCAR_interface.vasp')
    
