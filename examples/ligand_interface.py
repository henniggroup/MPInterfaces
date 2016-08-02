# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This script demonstrates the creation of ligand interfaces and 
preoptimization screening of possible interfaces.
"""

from six.moves import range

from pymatgen.core import Molecule, Structure
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.operations import SymmOp

from mpinterfaces.interface import Interface, Ligand
import numpy as np


def coloumb_configured_interface(iface, random=True,
                                 translations=None,
                                 rotations=None,
                                 samples=10, lowest=5, ecut=None):
    """
    Creates Ligand Slab interfaces of user specified translations 
    and rotations away from the initial guess of binding site 
    configuration, returns lowest energy structure according to 
    Coulomb model
    
    Args:
         Interface: Interface object: initial interface object 
         random: True for using Gaussian sampled random numbers for 
                 rotations and translations 
         translations: list of [x,y,z] translations to be performed
         rotation: list of [a,b,c] rotations to be performed w.r.t 
                   Ligand axis
         samples: number of interfaces to create
         lowest: number of structures to return according to order of 
                 minimum energies

    Returns: 
         list of lowest energy interface objects 
    """
    ifaces = []
    transform = []
    for i in range(samples):
        if random:
            x = np.random.normal()  # shift along x direction
            y = np.random.normal()  # shift along y direction
            z = np.random.normal()  # shift aling z direction
            a = SymmOp.from_axis_angle_and_translation(axis=[1, 0, 0], \
                                                       angle=np.random.normal(
                                                           0, 180))
            b = SymmOp.from_axis_angle_and_translation(axis=[0, 1, 0], \
                                                       angle=np.random.normal(
                                                           0, 180))
            c = SymmOp.from_axis_angle_and_translation(axis=[0, 0, 1], \
                                                       angle=np.random.normal(
                                                           0, 180))
            ligand = iface.ligand
            ligand.apply_operation(a)
            ligand.apply_operation(b)
            ligand.apply_operation(c)

        # check if created interface maintains the ligand adsorbed
        # over the surface 
        for j in iface.top_atoms:
            if not iface.cart_coords[j][2] + iface.displacement > \
                    min(ligand.cart_coords[:, 2]):
                transform.append(True)
        if all(transform):
            iface = Interface(iface.strt, hkl=iface.hkl,
                              min_thick=iface.min_thick,
                              min_vac=iface.min_vac,
                              supercell=iface.supercell,
                              surface_coverage=iface.surface_coverage,
                              ligand=iface.ligand, displacement=z,
                              adatom_on_lig=iface.adatom_on_lig,
                              adsorb_on_species=iface.adsorb_on_species,
                              primitive=False, from_ase=True,
                              x_shift=x, y_shift=y)
            iface.create_interface()
            energy = iface.calc_energy()
            iface.sort()
            if energy < ecut:
                ifaces.append((energy, iface))
                # ifaces.zip(energy, iface)
    return ifaces


if __name__ == '__main__':
    # PbS 100 surface with single hydrazine as ligand
    # sample pre-relaxed structure files for Bulk (strt) and Molecule
    strt = Structure.from_file("POSCAR_PbS")
    mol_struct = Structure.from_file("POSCAR_Hydrazine")
    mol = Molecule(mol_struct.species, mol_struct.cart_coords)
    hydrazine = Ligand([mol])

    # intital supercell, this wont be the final supercell if
    # surface coverage is specified
    supercell = [1, 1, 1]

    # miller index
    hkl = [1, 0, 0]

    # minimum slab thickness in Angstroms
    min_thick = 19

    # minimum vacuum thickness in Angstroms
    # mind: the ligand will be placed in this vacuum, so the
    # final effective vacuum space will be smaller than this
    min_vac = 12

    # surface coverage in the units of lig/ang^2
    # mind: exact coverage as provided cannot be guaranteed, the slab
    # will be constructed
    # with a coverage value thats close to the requested one
    # note: maximum supercell size possible is 10 x 10
    # note: 1 lig/nm^2 = 0.01 lig/ang^2
    surface_coverage = 0.01

    # atom on the slab surface on which the ligand will be attached,
    # no need to specify if the slab is made of only a single species
    adsorb_on_species = 'Pb'

    # atom on ligand that will be attached to the slab surface
    adatom_on_lig = 'N'

    # ligand displacement from the slab surface along the surface normal
    # i.e adatom_on_lig will be displced by this amount from the
    # adsorb_on_species atom
    # on the slab
    # in Angstrom
    displacement = 3.0

    # here we create the adsorbate slab Interface
    iface = Interface(strt, hkl=hkl, min_thick=min_thick,
                      min_vac=min_vac, supercell=supercell,
                      surface_coverage=surface_coverage,
                      ligand=hydrazine, displacement=displacement,
                      adatom_on_lig=adatom_on_lig,
                      adsorb_on_species=adsorb_on_species,
                      primitive=False, from_ase=True)
    iface.create_interface()
    iface.sort()
    energy = iface.calc_energy()
    iface.to('poscar', 'POSCAR_interface.vasp')
    interfaces = coloumb_configured_interface(iface, random=True,
                                              translations=None,
                                              rotations=None,
                                              samples=20, lowest=5,
                                              ecut=energy)
    for i, iface in enumerate(interfaces):
        print("Coloumb Energy")
        print(i, iface[0])
        iface[1].to('poscar', 'POSCAR_interface' + str(iface[0]) + '.vasp')
        iface_slab = iface[1].slab
        iface_slab.sort()
        # set selective dynamics flags as required
        true_site = [1, 1, 1]
        false_site = [0, 0, 0]
        sd_flag_iface = []
        sd_flag_slab = []
        # selective dynamics flags for the interface
        for i in iface.sites:
            sd_flag_iface.append(false_site)
        # selective dynamics flags for the bare slab
        for i in iface_slab.sites:
            sd_flag_slab.append(false_site)
            interface_poscar = Poscar(iface[1],
                                      selective_dynamics=sd_flag_iface)
            slab_poscar = Poscar(iface_slab,
                                 selective_dynamics=sd_flag_slab)
            # slab poscars without selective dynamics flag
            iface_slab.to('poscar', 'POSCAR_slab' + str(iface[0]) + '.vasp')
            # poscars with selective dynamics flag
            interface_poscar.write_file(
                "POSCAR_interface_with_sd" + str(iface[0]) + '.vasp')
            slab_poscar.write_file(
                "POSCAR_slab_with_sd" + str(iface[0]) + '.vasp')
