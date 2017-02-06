from __future__ import print_function, division, unicode_literals

import os

import yaml

from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.io.vasp.outputs import Vasprun
import twod_materials.utils as utl
from pymatgen.matproj.rest import MPRester

from monty.serialization import loadfn

import twod_materials


PACKAGE_PATH = twod_materials.__file__.replace('__init__.pyc', '')
PACKAGE_PATH = PACKAGE_PATH.replace('__init__.py', '')
PACKAGE_PATH = '/'.join(PACKAGE_PATH.split('/')[:-2])

try:
    config_vars = loadfn(os.path.join(os.path.expanduser('~'), 'config.yaml'))
except:
    print('WARNING: No config.yaml file was found. please configure the '\
    'config.yaml and put it in your home directory.')
    # Still set them for testing purposes.
    config_vars = loadfn(os.path.join(PACKAGE_PATH, 'config.yaml'))
if 'MP_API' in os.environ:  # Also for testing purposes.
    MPR = MPRester(os.environ['MP_API'])
else:
    MPR = MPRester(config_vars['mp_api'])
VASP = config_vars['normal_binary']
VASP_2D = config_vars['twod_binary']
if 'queue_system' in config_vars:
    QUEUE = config_vars['queue_system'].lower()
elif '/ufrc/' in os.getcwd():
    QUEUE = 'slurm'
elif '/scratch/' in os.getcwd():
    QUEUE = 'pbs'


class Calibrator():

    def __init__(self, incar_dict, potcar_dict, n_kpts_per_atom=500,
                 ncores=1, nprocs=16, pmem='600mb', walltime='6:00:00',
                 binary='vasp'):
        """
        Args:
            incar_dict (dict): dictionary of all input parameters
                used in the given framework.
            potcar_dict (dict): dictionary of all species to be
                calibrated and the potcar hashes used in the
                given framework, e.g. {'Mo': 'pv', 'S': ''}.
            n_kpts_per_atom (int): Create kpoints at specified
                density per atom. Defaults to 500.
            n_cores, n_procs, pmem, walltime, binary: runjob
                parameters. Defaults established for a regular
                sized job on hipergator.
        """

        self._incar_dict = incar_dict
        self._n_kpts_per_atom = n_kpts_per_atom
        self._potcar_dict = potcar_dict
        self._ncores = ncores
        self._nprocs = nprocs
        self._pmem = pmem
        self._walltime = walltime
        self._binary = binary
        self._config = loadfn('/home/mashton/cal_config.yaml')

    def prepare(self, submit=False):
        """
        Set up calculation directories to calibrate
        the ion corrections to match a specified framework of INCAR
        parameters, k-points, and potcar hashes.

        Args:
            submit (bool): whether or not to submit each job
                after preparing it.
        """

        for elt in self._potcar_dict:

            # Set up reference directory for the pure element.
            if not os.path.isdir(elt):
                os.mkdir(elt)
            os.chdir(elt)

            # Poscar
            s = MPR.get_structure_by_material_id(
                self._config['Mpids'][elt]['self']
                )
            s.to('POSCAR', 'POSCAR')
            plines = open('POSCAR').readlines()
            elements = plines[5].split()

            # Kpoints
            kp = Kpoints.automatic_density(s, self._n_kpts_per_atom)
            kp.write_file('KPOINTS')

            # Incar
            incar = Incar.from_dict(self._incar_dict)
            incar.write_file('INCAR')

            # Potcar
            utl.write_potcar(types=[self._potcar_dict[el] for el in elements])

            # Runjob

            if QUEUE == 'pbs':
                utl.write_pbs_runjob('{}_cal'.format(elt), self._ncores,
                                     self._nprocs, self._pmem, self._walltime,
                                     self._binary)
                submission_command = 'qsub runjob'

            elif QUEUE == 'slurm':
                utl.write_slurm_runjob('{}_cal'.format(elt), self._nprocs,
                                       self._pmem, self._walltime,
                                       self._binary)
                submission_command = 'sbatch runjob'

            if submit:
                os.system(submission_command)

            # Set up reference oxide compound subdirectory.
            if elt not in ['O', 'S', 'F', 'Cl', 'Br', 'I']:
                if not os.path.isdir('ref'):
                    os.mkdir('ref')
                os.chdir('ref')

                # Poscar
                s = MPR.get_structure_by_material_id(
                    self._config['Mpids'][elt]['ref']
                    )
                s.to('POSCAR', 'POSCAR')
                plines = open('POSCAR').readlines()
                elements = plines[5].split()

                # Kpoints
                kp = Kpoints.automatic_density(s, self._n_kpts_per_atom)
                kp.write_file('KPOINTS')

                # Incar
                incar = Incar.from_dict(self._incar_dict)
                incar.write_file('INCAR')

                # Potcar
                utl.write_potcar(
                    types=[self._potcar_dict[el] for el in elements])

                # Runjob
                if QUEUE == 'slurm':
                    utl.write_pbs_runjob('{}_cal'.format(elt), self._ncores,
                                         self._nprocs, self._pmem,
                                         self._walltime, self._binary)
                    submission_command = 'qsub runjob'

                elif QUEUE == 'pbs':
                    utl.write_slurm_runjob('{}_cal'.format(elt), self._nprocs,
                                           self._pmem, self._walltime,
                                           self._binary)
                    submission_command = 'sbatch runjob'

                if submit:
                    os.system(submission_command)

                os.chdir('../')
            os.chdir('../')

    def get_corrections(self, parent_dir=os.getcwd(), write_yaml=False,
                        oxide_corr=0.708):
        """
        Pulls the corrections to be added for each element.

        Args:
            parent_dir (str): path to parent directory containing
                subdirectories created by prepare(). Defaults to cwd.
            write_yaml (bool): whether or not to write the
                corrections to ion_corrections.yaml and the mu0
                values to end_members.yaml.
            oxide_corr (float): additional correction added for oxygen
                to get water's formation energy right.

        Returns:
            dict. elements as keys and their corrections as values,
                in eV per atom, e.g. {'Mo': 0.135, 'S': -0.664}.
        """

        mu0 = dict()
        corrections = dict()

        os.chdir(parent_dir)

        special_cases = ['O', 'S', 'F', 'Cl', 'Br', 'I']

        elts = [elt for elt in self._potcar_dict if elt not in special_cases]

        # Add entropic correction for special elements (S * 298K)
        specials = [elt for elt in self._potcar_dict if elt in special_cases]
        for elt in specials:
            os.chdir(elt)
            vasprun = Vasprun('vasprun.xml')
            composition = vasprun.final_structure.composition
            n_formula_units = composition.get_integer_formula_and_factor()[1]

            mu0[elt] = (
                round(vasprun.final_energy / n_formula_units
                      + self._config['OtherCorrections'][elt], 3)
                )
            os.chdir(parent_dir)

        # Oxide correction from Materials Project
        mu0['O'] += oxide_corr

        for elt in elts:
            os.chdir(elt)
            vasprun = Vasprun('vasprun.xml')
            composition = vasprun.final_structure.composition
            n_formula_units = composition.get_integer_formula_and_factor()[1]

            mu0[elt] = round(vasprun.final_energy / n_formula_units, 3)

            # Nitrogen needs both kinds of corrections
            if elt == 'N':
                mu0[elt] -= 0.296

            os.chdir(parent_dir)

        for elt in elts:

            os.chdir('{}/ref'.format(elt))
            vasprun = Vasprun('vasprun.xml')
            composition = vasprun.final_structure.composition
            n_formula_units = composition.get_integer_formula_and_factor()[1]

            fH_exp = self._config['Experimental_fH'][elt]
            try:
                fH_dft = vasprun.final_energy / n_formula_units

                plines = open('POSCAR').readlines()
                elements = plines[5].split()
                stoichiometries = plines[6].split()
                comp_as_dict = {}
                for element in elements:
                    comp_as_dict[element] = 0
                for i, element in enumerate(elements):
                    comp_as_dict[element] += int(stoichiometries[i])

                n_elt_per_fu = (
                    int(comp_as_dict[elt]) / n_formula_units
                    )
                for el in comp_as_dict:
                    fH_dft -= (
                        mu0[el] * int(comp_as_dict[el])
                        / n_formula_units
                        )

                corrections[elt] = round((fH_dft - fH_exp) / n_elt_per_fu, 3)

            except UnboundLocalError:
                corrections[elt] = 'Not finished'

            os.chdir(parent_dir)

        if write_yaml:
            with open('ion_corrections.yaml', 'w') as icy:
                icy.write(yaml.dump(corrections, default_flow_style=False))
            with open('end_members.yaml', 'w') as emy:
                emy.write(yaml.dump(mu0, default_flow_style=False))

        return corrections
