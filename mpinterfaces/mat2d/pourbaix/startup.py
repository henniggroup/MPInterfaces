from __future__ import print_function, division, unicode_literals

import os
import yaml

from monty.serialization import loadfn

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.outputs import Vasprun

from mpinterfaces import MPR
import mpinterfaces.utils as utl
from mpinterfaces.mat2d.stability import relax

__author__ = "Michael Ashton"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Michael Ashton"
__email__ = "ashtonmv@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

PACKAGE_PATH = os.path.dirname(__file__)
REFERENCE_MPIDS = loadfn(os.path.join(PACKAGE_PATH, 'reference_mpids.yaml'))
EXPERIMENTAL_OXIDE_DATA = loadfn(os.path.join(PACKAGE_PATH,
                                              'experimental_oxide_data.yaml'))
GAS_CORRECTIONS = loadfn(os.path.join(PACKAGE_PATH, 'gas_corrections.yaml'))


def get_experimental_formation_energies():
    """
    Read in the raw enthalpy and entropy energy data from
    Kubaschewski in experimental_oxide_data.yaml
    and interpret it into actual formation energies. This extra
    step is written out mostly just to make the methodology clear
    and reproducible.

    Returns:

    """
    data = EXPERIMENTAL_OXIDE_DATA
    oxygen_entropy = 38.48  # Entropy in atomic O, in cal/mol.degree
    formation_energies = {}
    for compound in data:
        composition = Composition(compound)
        element = [e for e in composition if e.symbol != 'O'][0]
        delta_H = data[compound]['delta_H']
        delta_S = (data[compound]['S_cmpd']
                   - (data[compound]['S_elt']*composition[element]
                      + oxygen_entropy*composition['O'])) * 298 / 1000
        # Convert kcal/mole to eV/formula unit
        formation_energies[element.symbol] = (delta_H - delta_S) / 22.06035

    return formation_energies


def relax_references(potcar_types, incar_dict, submit=True,
                     force_overwrite=False):
    """
    Set up calculation directories to calibrate
    the ion corrections to match a specified framework of INCAR
    parameters and potcar hashes.

    Args:
        potcar_types (list): list of all elements to calibrate,
            containing specifications for the kind of potential
            desired for each element, e.g. ['Na_pv', 'O']. If
            oxygen is not explicitly included in the list, 'O'
            is used.

        incar_dict (dict): a dictionary of input parameters
            and their values, e.g. {'ISMEAR': -5, 'NSW': 10}

        submit (bool): whether or not to submit each job
            after preparing it.

        force_overwrite (bool): Whether or not to overwrite files
            if an already converged vasprun.xml exists in the
            directory.
    """

    for element in potcar_types:
        if element.split('_')[0] == 'O':
            oxygen_potcar = element
            break
    else:
        oxygen_potcar = 'O'
        potcar_types.append('O')

    for element in potcar_types:
        elt = element.split('_')[0]

        # First, set up a relaxation for the pure element.
        if not os.path.isdir(elt):
            os.mkdir(elt)
        os.chdir(elt)
        s = MPR.get_structure_by_material_id(REFERENCE_MPIDS[elt]['element'])
        s.to('POSCAR', 'POSCAR')
        relax(dim=3, submit=submit, force_overwrite=force_overwrite)
        incar_dict.update({"MAGMOM": utl.get_magmom_string(s)})
        Incar.from_dict(incar_dict).write_file("INCAR")
        utl.write_potcar(types=[element])

        # Then set up a relaxation for its reference oxide.
        if elt not in ['O', 'S', 'F', 'Cl', 'Br', 'I', 'H']:
            if not os.path.isdir('oxide'):
                os.mkdir('oxide')
            os.chdir('oxide')
            s = MPR.get_structure_by_material_id(REFERENCE_MPIDS[elt]['oxide'])
            s.get_sorted_structure().to('POSCAR', 'POSCAR')
            relax(dim=3, submit=submit, force_overwrite=force_overwrite)
            incar_dict.update({"MAGMOM": utl.get_magmom_string(s)})
            Incar.from_dict(incar_dict).write_file("INCAR")
            utl.write_potcar(types=[element, oxygen_potcar])

            os.chdir('../')
        os.chdir('../')


def get_corrections(write_yaml=False):
    """
    Calculates and collects the corrections to be added for
    each reference element directory in the current working
    directory.

    Args:
        write_yaml (bool): whether or not to write the
            corrections to ion_corrections.yaml and the mu0
            values to end_members.yaml.

    Returns:
        dict: elements as keys and their corrections as values,
            in eV per atom, e.g. {'Mo': 0.135, 'S': -0.664}.
    """

    experimental_formation_energies = get_experimental_formation_energies()
    mu0, corrections = {}, {}
    special_cases = ['O', 'S', 'F', 'Cl', 'Br', 'I', 'H']

    elts = [elt for elt in os.listdir(os.getcwd())
            if os.path.isdir(elt) and elt not in special_cases]
    special_elts = [elt for elt in os.listdir(os.getcwd())
                    if os.path.isdir(elt) and elt in special_cases]

    # Add entropic correction for special elements (S * 298K)
    for elt in special_elts:
        os.chdir(elt)
        try:
            vasprun = Vasprun('vasprun.xml')
            composition = vasprun.final_structure.composition
            formula_and_factor = composition.get_integer_formula_and_factor()
            n_formula_units = composition.get_integer_formula_and_factor()[1]
            if '2' in formula_and_factor[0]:
                n_formula_units *= 2

            mu0[elt] = (round(vasprun.final_energy / n_formula_units
                              + GAS_CORRECTIONS['entropy'][elt], 3))
        except Exception as e:
            mu0[elt] = 'Element not finished'
        os.chdir('../')

    # Oxide correction based on L. Wang, T. Maxisch, and G. Ceder,
    # Phys. Rev. B 73, 195107 (2006). This correction is to get
    # solid oxide formation energies right.
    #mu0['O'] += GAS_CORRECTIONS['oxide']

    for elt in elts:
        EF_exp = experimental_formation_energies[elt]
        os.chdir(elt)
        try:
            vasprun = Vasprun('vasprun.xml')
            composition = vasprun.final_structure.composition
            mu0[elt] = round(vasprun.final_energy/composition[Element(elt)], 3)

            # Nitrogen needs an entropic gas phase correction too.
            if elt == 'N':
                mu0[elt] -= GAS_CORRECTIONS['entropy']['N']

        except Exception as _:
            corrections[elt] = 'Element not finished'

        os.chdir('oxide')
        try:
            vasprun = Vasprun('vasprun.xml')
            composition = vasprun.final_structure.composition
            n_formula_units = composition.get_integer_formula_and_factor()[1]

            EF_dft = (vasprun.final_energy
                      - mu0[elt]*composition[Element(elt)]
                      - mu0['O']*composition[Element('O')]) / n_formula_units

            corrections[elt] = round((EF_dft - EF_exp)
                                     / (composition[Element(elt)]/n_formula_units), 3)

        except Exception as e:
            print(e)
            # The relaxation didn't finish.
            if elt in corrections:
                corrections[elt] += ' and oxide not finished'
            else:
                corrections[elt] = 'Oxide not finished'

        os.chdir('../../')
    if write_yaml:
        with open('ion_corrections.yaml', 'w') as yam:
            yam.write('# Difference in formation energy between\n')
            yam.write('# DFT framework and experimental data\n')
            yam.write('# from Materials Thermochemistry\n')
            yam.write('# (Kubaschewski, Alcock)\n')
            yam.write('# All values are given in eV per atom.\n\n')
            yam.write(yaml.dump(corrections, default_flow_style=False))
        with open('chemical_potentials.yaml', 'w') as yam:
            yam.write('# chemical potentials of elemental phases to calculate')
            yam.write('\n')
            yam.write('# formation energies of all compounds. In eV/atom.\n\n')
            yam.write(yaml.dump(mu0, default_flow_style=False))

    return corrections
