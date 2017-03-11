# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Wulff construction to create the nanoparticle
"""

from six.moves import range

import itertools
from math import gcd
from functools import reduce

import numpy as np

from pymatgen.core.structure import Structure, Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord_utils import in_coord_list

from mpinterfaces import get_struct_from_mp

from mpinterfaces.default_logger import get_default_logger

logger = get_default_logger(__name__)


class Nanoparticle(Molecule):
    """
    Construct nanoparticle using wulff construction
    """

    def __init__(self, structure, rmax=15, hkl_family=((1, 0, 0), (1, 1, 1)),
                 surface_energies=(28, 25)):
        self.structure = structure
        self.rmax = rmax
        self.hkl_family = list(hkl_family)
        self.surface_energies = list(surface_energies)
        spherical_neighbors = self.structure.get_sites_in_sphere(
            [0.0, 0.0, 0.0], self.rmax)
        recp_lattice = self.structure.lattice.reciprocal_lattice_crystallographic
        self.recp_lattice = recp_lattice.scale(1)
        self.set_miller_family()
        Molecule.__init__(self, [sn[0].species_and_occu
                                 for sn in spherical_neighbors],
                          [sn[0].coords for sn in spherical_neighbors],
                          charge=0)

    def set_miller_family(self):
        """
        get all miller indices for the given maximum index
        get the list of indices that correspond to the given family
        of indices
        """
        recp_structure = Structure(self.recp_lattice, ["H"], [[0, 0, 0]])
        analyzer = SpacegroupAnalyzer(recp_structure, symprec=0.001)
        symm_ops = analyzer.get_symmetry_operations()

        max_index = max(max(m) for m in self.hkl_family)
        r = list(range(-max_index, max_index + 1))
        r.reverse()
        miller_indices = []
        self.all_equiv_millers = []
        self.all_surface_energies = []
        for miller in itertools.product(r, r, r):
            if any([i != 0 for i in miller]):
                d = abs(reduce(gcd, miller))
                miller_index = tuple([int(i / d) for i in miller])
                for op in symm_ops:
                    for i, u_miller in enumerate(self.hkl_family):
                        if in_coord_list(u_miller, op.operate(miller_index)):
                            self.all_equiv_millers.append(miller_index)
                            self.all_surface_energies.append(
                                self.surface_energies[i])

    def get_normals(self):
        """
        get the normal to the plane (h,k,l)
        """
        normals = []
        for hkl in self.all_equiv_millers:
            normal = self.recp_lattice.matrix[0, :] * hkl[0] + \
                self.recp_lattice.matrix[1, :] * hkl[1] + \
                self.recp_lattice.matrix[2, :] * hkl[2]
            normals.append(normal / np.linalg.norm(normal))
        return normals

    def get_centered_molecule(self):
        center = self.center_of_mass
        new_coords = np.array(self.cart_coords) - center
        return Molecule(self.species_and_occu, new_coords,
                        charge=self._charge,
                        spin_multiplicity=self._spin_multiplicity,
                        site_properties=self.site_properties)

    def create(self):
        """
        creates the nanoparticle by chopping of the corners normal to the
        specified surfaces.
        the distance to the surface from the center of the particel =
        normalized surface energy * max radius
        """
        mol = self.get_centered_molecule()
        normalized_surface_energies = \
            np.array(self.all_surface_energies) / float(
                max(self.all_surface_energies))
        surface_normals = self.get_normals()
        remove_sites = []
        for i, site in enumerate(mol):
            for j, normal in enumerate(surface_normals):
                n = np.array(normal)
                n = n / np.linalg.norm(n)
                if np.dot(site.coords, n) + self.rmax * \
                        normalized_surface_energies[j] <= 0:
                    remove_sites.append(i)
                    break
        self.remove_sites(remove_sites)
        # new_sites = [site for k, site in enumerate(mol) if k not in remove_sites]
        # return Molecule.from_sites(new_sites)


if __name__ == '__main__':
    # nanopartcle settings
    # max radius in angstroms
    rmax = 15
    # surface families to be chopped off
    surface_families = [(1, 0, 0), (1, 1, 1)]
    # could be in any units, will be normalized
    surface_energies = [28, 25]

    # caution: set the structure wrt which the the miller indices are specified
    # use your own API key
    structure = get_struct_from_mp('PbS')
    # primitve ---> conventional cell
    sa = SpacegroupAnalyzer(structure)
    structure_conventional = sa.get_conventional_standard_structure()

    nanoparticle = Nanoparticle(structure_conventional, rmax=rmax,
                                hkl_family=surface_families,
                                surface_energies=surface_energies)
    nanoparticle.create()
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
