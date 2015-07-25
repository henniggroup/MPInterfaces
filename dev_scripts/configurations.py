from __future__ import division, unicode_literals, print_function

"""
Compute all configurations of 2D materials on substrates by identifying
the Wyckoff positions of the atoms.
To be used after determining the super-lattices from the hetero_inteface.py
script
"""

import sys 

import numpy as np

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure


def get_uniq_layercoords(struct, nlayers, top=True):
    """
    returns the coordinates of unique sites in the top or bottom
    nlayers of the given structure.

    Args:
        struct: input structure
        nlayers: number of layers
        top: top or bottom layers, default is top layer

    Return:
        numpy array of unique coordinates
    """
    coords = np.array([site.coords for site in struct])
    z = coords[:,2]
    zu, zuind, = np.unique(z, return_index=True)
    if top:
        z_nthlayer = z[zuind[-nlayers]]
        zfilter = (z >= z_nthlayer)
    else:
        z_nthlayer = z[zuind[nlayers - 1]]
        zfilter = (z <= z_nthlayer)
    # site indices in the layers        
    indices_layers = np.argwhere(zfilter).ravel()
    sa = SpacegroupAnalyzer(struct)
    symm_data = sa.get_symmetry_dataset()
    # equivalency mapping for the structure
    # i'th site in the struct equivalent to eq_struct[i]'th site
    eq_struct = symm_data["equivalent_atoms"]
    # equivalency mapping for the layers
    eq_layers = eq_struct[indices_layers]
    print(eq_layers)
    # site indices of unique atoms in the layers
    __, ueq_layers_indices = np.unique(eq_layers, return_index=True)
    print(ueq_layers_indices) 
    indices_uniq = indices_layers[ueq_layers_indices]
    print(indices_uniq)
    # coordinates of the unique atoms in the layers
    print(coords[indices_uniq])
    return coords[indices_uniq]
    

def generate_all_configs(mat2d, substrate,
                         nlayers_2d = 2, nlayers_substrate = 2 ,
                         seperation = 5 ):
    """
    For the given lattice matched 2D material and substrate structures,
    this functions computes all unique sites in the interface layers
    and subsequently generates all possible unique 2d/substrate
    interfaces and writes the corresponding poscar files

    Args:
        mat2d: Lattice and symmetry-matched 2D material structure
        substrate: Lattice and symmetry-matched 2D substrate structure
        nlayers_substrate: number of substrate layers
        nlayers_2d: number of 2d material layers
        seperation: seperation between the substrate and the 2d
                    material
        
    Returns:
        None
        
    TODO: give additional random placement of 2D material on substrate
    """
    # immediate exit if no structures
    if not(mat2d and substrate):
        print("no structures. aborting ...")
        sys.exit()    
    # unique site coordinates in the substrate top layers
    coords_uniq_sub = get_uniq_layercoords(substrate,
                                           nlayers_substrate,
                                           top=True)
    # unique site coordinates in the 2D material bottom layers    
    coords_uniq_2d = get_uniq_layercoords(mat2d,
                                          nlayers_2d,
                                          top=False)
    substrate_top_z = np.max(np.array([site.coords
                                       for site in substrate])[:,2])
    mat_2d_bottom = np.min(np.array([site.coords
                                     for site in mat2d])[:,2])
    # shift normal to the surface by 'seperation'
    surface_normal = substrate.lattice.matrix[2,:]
    origin = np.array([0,0,substrate_top_z])
    shift_normal = surface_normal/np.linalg.norm(surface_normal) * seperation

    #generate all possible interfaces, one for each combination of
    # unique substrate and unique 2d materials site in the layers .i.e
    # an interface structure for each parallel shift
    #interface = 2D material + substrate
    for i, coord_i in enumerate(coords_uniq_sub):
        for j, coord_j in enumerate(coords_uniq_2d):
            interface = substrate.copy()
            shift_parallel = coord_i - coord_j
            shift_parallel[2] = 0
            shift_net = shift_normal - shift_parallel
            for site in mat2d:
                new_coords = site.coords
                new_coords[2] = site.coords[2] - mat_2d_bottom
                new_coords = new_coords + origin + shift_net
                interface.append(site.specie, new_coords,
                                 coords_are_cartesian=True)
            interface.to(fmt='poscar',
                         filename='POSCAR_final_'+str(i)+'_'+str(j)+'.vasp')

            
# test
if __name__ == '__main__':            
    mat2d = Structure.from_file('POSCAR_2D')
    substrate = Structure.from_file('POSCAR_substrate')
    seperation = 5
    nlayers_2d = 2
    nlayers_substrate = 2
    generate_all_configs(mat2d, substrate,
                       nlayers_2d, nlayers_substrate,
                       seperation )
