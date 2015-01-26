"""
Wulff construction to create the nanoparticle
"""

import sys
import itertools
from fractions import gcd
from functools import reduce
from collections import defaultdict

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.surface import get_symmetrically_distinct_miller_indices
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord_utils import in_coord_list

from mpinterfaces import get_struct_from_mp

def get_normals(recp_lattice, hkls):
    """
    get the normal to the plane (h,k,l)
    """
    normals = []
    for hkl in hkls:
        normal = recp_lattice.matrix[0,:]*hkl[0] + \
          recp_lattice.matrix[1,:]*hkl[1] + \
          recp_lattice.matrix[2,:]*hkl[2] 
        normals.append(normal/np.linalg.norm(normal))
    return normals

def create_nanoparticle(recp_lattice, structure, hkls, surface_energies, rmax):
    """
    creates the nanoparticle by chopping of the corners normal to the
    specified surfaces.
    the distance to the surface from the center of the particel =
    normalized surface energy * max radius
    """
    mol = Molecule( [site.species_and_occu
                              for site in structure.sites],
                              structure.cart_coords)
    mol = mol.get_centered_molecule()    
    surface_energies = np.array(surface_energies)/ float(max(surface_energies))
    surface_normals = get_normals(recp_lattice, hkls) #hkls
    remove_sites = []
    for i, site in enumerate(mol):
        for j, normal in enumerate(surface_normals):
            normal = surface_normals[j]
            surface_energy = surface_energies[j]
            n = np.array(normal)
            n = n/np.linalg.norm(n)
            if np.dot(site.coords, n) + rmax*surface_energy <= 0 :
                remove_sites.append(i)
                break
    new_sites = [site for k, site in enumerate(mol)
                 if k not in remove_sites]
    return Molecule.from_sites(new_sites)

            
if __name__ == '__main__':
    #nanopartcle settings
    rmax = 15 #max radius in angstroms
    surfaces = [(1,0,0), (1,1,1)] #surface families to be chopped off
    surface_energies = [28,25] #could be in any units, will be normalized 

    #caution: set the structure wrt which the the miller indices are specified
    #use your own key
    structure = get_struct_from_mp('PbS', MAPI_KEY="dwvz2XCFUEI9fJiR")
    #primitve ---> conventional cell
    sa = SpacegroupAnalyzer(structure)
    structure_conventional = sa.get_conventional_standard_structure()

    #supercell from which the nanoparticle is carved
    ncella = int(np.ceil(2*rmax/structure_conventional.lattice.a))
    ncellb = int(np.ceil(2*rmax/structure_conventional.lattice.b))
    ncellc = int(np.ceil(2*rmax/structure_conventional.lattice.c))        
    structure_conventional.make_supercell([ncella, ncellb, ncellc])
    structure_conventional.to(fmt='poscar', filename='nano_supercell.vasp')    
    
    #uniq_millers = get_symmetrically_distinct_miller_indices(structure,1)
    recp_lattice = structure_conventional.lattice.reciprocal_lattice_crystallographic
    recp_lattice = recp_lattice.scale(1)
    recp = Structure(recp_lattice, ["H"], [[0, 0, 0]])
    analyzer = SpacegroupAnalyzer(recp, symprec=0.001)
    symm_ops = analyzer.get_symmetry_operations()

    #get all miller indices for the given maximum index
    #get the list of indices that correspond to the given family
    #of indices    
    max_index = max( max(m) for m in surfaces)
    r = list(range(-max_index, max_index + 1))
    r.reverse()
    miller_indices = []
    equiv_millers = []
    s_energies = []    
    for miller in itertools.product(r, r, r):
        if any([i != 0 for i in miller]):
            d = abs(reduce(gcd, miller))
            miller_index = tuple([int(i / d) for i in miller])
            for op in symm_ops:
                for i, u_miller in enumerate(surfaces):
                    if in_coord_list(u_miller, op.operate(miller_index)):
                        equiv_millers.append(miller_index)
                        s_energies.append(surface_energies[i])


    nanoparticle = create_nanoparticle(recp_lattice, structure_conventional,
                                       equiv_millers, s_energies, rmax)
    nanoparticle.to(fmt='xyz', filename='nanoparticle.xyz')

    """
    Wulff construction using the ASE package
    works only for cubic systems and doesn't support multiatom basis
    
    from ase.cluster import wulff_construction
    from pymatgen.io.aseio import AseAtomsAdaptor

    symbol = 'Pt'
    surfaces = [ (1,0,0), (1,1,1) ]
    surface_energies = [1, 1]
    size = 200 #number of atoms
    structure = "fcc"
    latticeconstant = 5.0
    atoms = wulff_construction(symbol, surfaces, surface_energies, size, structure,
    rounding='closest', latticeconstant=latticeconstant,
    debug=False, maxiter=100)
    #convert to pymatgen structure
    pgen_structure = AseAtomsAdaptor().get_structure(atoms)
    pgen_structure.to(fmt='poscar', filename='POSCAR_pt_nano.vasp')
    
    """
