from __future__ import print_function, division, unicode_literals

import os

import yaml

from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.matproj.rest import MPRester

from monty.serialization import loadfn

from mpinterfaces import MY_CONFIG

from twod_materials import MPR, VASP, VASP_2D, POTENTIAL_PATH, USR, VDW_KERNEL,\
    QUEUE
import twod_materials
from twod_materials.stability.startup import relax
import twod_materials.utils as utl


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')

REFERENCES = loadfn(os.path.join(
    PACKAGE_PATH, 'pourbaix/reference_oxides.yaml')
)


def relax_references(potcar_types, incar_dict, submit=True):
    """
    Set up calculation directories to calibrate
    the ion corrections to match a specified framework of INCAR
    parameters and potcar hashes.

    Args:
        potcar_types (list): list of all elements to calibrate,
            containing specifications for the kind of potential
            desired for each element, e.g. ['Na_pv', 'O_s']. If
            oxygen is not explicitly included in the list, 'O_s'
            is used.

        incar_dict (dict): a dictionary of input parameters
            and their values, e.g. {'ISMEAR': -5, 'NSW': 10}

        submit (bool): whether or not to submit each job
            after preparing it.
    """

    for element in potcar_types:
        if element.split('_')[0] == 'O':
            oxygen_potcar = element
            break
    else:
        oxygen_potcar = 'O_s'

    for element in potcar_types:
        elt = element.split('_')[0]
        # First, set up a relaxation for the pure element.
        if not os.path.isdir(elt):
            os.mkdir(elt)
        os.chdir(elt)
        s = MPR.get_structure_by_material_id(
            REFERENCES['Mpids'][elt]['self']
            )
        s.to('POSCAR', 'POSCAR')
        relax(dim=3, incar_dict=incar_dict, submit=submit)
        utl.write_potcar(types=[element])

        # Then set up a relaxation for its reference oxide.
        if elt not in ['O', 'S', 'F', 'Cl', 'Br', 'I']:
            if not os.path.isdir('ref'):
                os.mkdir('ref')
            os.chdir('ref')
            s = MPR.get_structure_by_material_id(
                REFERENCES['Mpids'][elt]['ref']
                )
            s.to('POSCAR', 'POSCAR')
            relax(dim=3, incar_dict=incar_dict, submit=submit)
            utl.write_potcar(types=[element, oxygen_potcar])

            os.chdir('../')
        os.chdir('../')


def get_corrections(write_yaml=False, oxide_corr=0.708):
    """
    Calculates and collects the corrections to be added for
    each reference element directory in the current working
    directory.

    Args:
        write_yaml (bool): whether or not to write the
            corrections to ion_corrections.yaml and the mu0
            values to end_members.yaml.

        oxide_corr (float): additional correction added for oxygen
            to get water's formation energy right.

    Returns:
        dict. elements as keys and their corrections as values,
            in eV per atom, e.g. {'Mo': 0.135, 'S': -0.664}.
    """

    mu0, corrections = {}, {}

    special_cases = ['O', 'S', 'F', 'Cl', 'Br', 'I']

    elts = [elt for elt in os.listdir(os.getcwd()) if os.path.isdir(elt)
            and elt not in special_cases]
    special_elts = [elt for elt in os.listdir(os.getcwd()) if os.path.isdir(elt)
            and elt in special_cases]

    # Add entropic correction for special elements (S * 298K)
    for elt in special_elts:
        os.chdir(elt)
        vasprun = Vasprun('vasprun.xml')
        composition = vasprun.final_structure.composition
        n_formula_units = composition.get_integer_formula_and_factor()[1]

        mu0[elt] = (
            round(vasprun.final_energy / n_formula_units
                  + REFERENCES['OtherCorrections'][elt], 3)
            )
        os.chdir(parent_dir)

    # Oxide correction from Materials Project
    mu0['O'] += oxide_corr

    for elt in elts:
        fH_exp = REFERENCES['Experimental_fH'][elt]

        os.chdir(elt)
        try:
            vasprun = Vasprun('vasprun.xml')
            composition = vasprun.final_structure.composition
            n_formula_units = composition.get_integer_formula_and_factor()[1]

            mu0[elt] = round(vasprun.final_energy / n_formula_units, 3)

            # Nitrogen needs an entropic gas phase correction too.
            if elt == 'N':
                mu0[elt] -= 0.296
        except Exception as e:
            corrections[elt] = 'Element not finished'

        os.chdir('ref')
        try:
            vasprun = Vasprun('vasprun.xml')
            composition = vasprun.final_structure.composition
            n_formula_units = composition.get_integer_formula_and_factor()[1]

            fH_dft = vasprun.final_energy / n_formula_units
            n_elt_per_formula_unit = composition[Element(elt)]
            corrections[elt] = round((fH_dft - fH_exp) / n_elt_per_fu, 3)

        except UnboundLocalError:
            # The relaxation didn't finish.
            if elt in corrections:
                corrections[elt] += 'and oxide not finished'
            else:
                corrections[elt] = 'Oxide not finished'

        os.chdir('../../')

    if write_yaml:
        with open('ion_corrections.yaml', 'w') as icy:
            icy.write(yaml.dump(corrections, default_flow_style=False))
        with open('end_members.yaml', 'w') as emy:
            emy.write(yaml.dump(mu0, default_flow_style=False))

    return corrections
