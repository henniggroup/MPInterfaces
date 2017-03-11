# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Defines measurement jobs
"""

from six.moves import zip

import sys
import shutil
import os
import json
import itertools

from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.inputs import Potcar, Kpoints

from mpinterfaces.calibrate import CalibrateMolecule
from mpinterfaces.calibrate import CalibrateSlab
from mpinterfaces.calibrate import CalibrateInterface
from mpinterfaces.interface import Interface
from mpinterfaces.default_logger import get_default_logger

__author__ = "Kiran Mathew, Joshua J. Gabriel"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

logger = get_default_logger(__name__)


class Measurement(object):
    """
    Takes in calibrate objects and use that to perform various
    measurement calculations such as solvation, ligand binding energy
    etc. The default behaviour is to setup and run static calculations
    for all the given calibrate jobs.

    Serves as Base Class. Override this class for custom measurements.

    Args:
        cal_objs: List of Calibration Object Names
        parent_job_dir: Directory in which measuremnt is setup
        job_dir: Path name to directory for running the Measurement
                 modules
    """

    def __init__(self, cal_objs, parent_job_dir='.', job_dir='./Measurement'):
        self.jobs = []
        self.handlers = []
        self.calmol = []
        self.calslab = []
        self.calbulk = []
        self.cal_objs = cal_objs
        self.job_dir = job_dir
        for obj in cal_objs:
            obj.old_jobs = obj.jobs
            obj.jobs = []
            obj.old_job_dir_list = obj.job_dir_list
            obj.job_dir_list = []

    def setup(self):
        """
        setup static jobs for all the calibrate objects
        copies CONTCAR to POSCAR
        sets NSW = 0
        """
        for cal in self.cal_objs:
            for i, jdir in enumerate(cal.old_job_dir_list):
                job_dir = self.job_dir + os.sep \
                    + jdir.replace(os.sep, '_').replace('.', '_') \
                    + os.sep + 'STATIC'
                logger.info('setting up job in {}'.format(job_dir))
                cal.incar = Incar.from_file(jdir + os.sep + 'INCAR')
                cal.incar['EDIFF'] = '1E-6'
                cal.incar['NSW'] = 0
                cal.potcar = Potcar.from_file(jdir + os.sep + 'POTCAR')
                cal.kpoints = Kpoints.from_file(jdir + os.sep + 'KPOINTS')
                contcar_file = jdir + os.sep + 'CONTCAR'
                if os.path.isfile(contcar_file):
                    logger.info('setting poscar file from {}'
                                .format(contcar_file))
                    cal.poscar = Poscar.from_file(contcar_file)
                    cal.add_job(job_dir=job_dir)
                else:
                    logger.critical("""CONTCAR doesnt exist.
                    Setting up job using input set in the old
                    calibration directory""")
                    cal.poscar = Poscar.from_file(jdir + os.sep + 'POSCAR')
                    cal.add_job(job_dir=job_dir)

    def run(self, job_cmd=None):
        """ run jobs """
        for cal in self.cal_objs:
            cal.run(job_cmd=job_cmd)

    def get_energy(self, cal):
        """
         measures the energy of a single cal object
         a single cal object can have multiple calculations
         returns energies lists
        """

    pass

    #     energies = []
    #     for job_dir in cal.job_dir_list:
    #         drone = MPINTVaspDrone(inc_structure=True,
    #                                inc_incar_n_kpoints=False)
    #         bg = BorgQueen(drone)
    #         # bg.parallel_assimilate(rootpath)
    #         bg.serial_assimilate(job_dir)
    #         allentries = bg.get_data()
    #         for e in allentries:
    #             if e:
    #                 energies.append(e.energy)
    #                 logger.debug('energy from directory {0} : {1}'
    #                              .format(job_dir, e.energy))
    #     return energies

    def make_measurements(self):
        """
        gets the energies and processes it
        override this for custom measurements
        """
        energies = []
        for cal in self.cal_objs:
            energies.append(self.get_energy(cal))


class MeasurementSolvation(Measurement):
    """
    Solvation with poisson-boltzmann(test verison)

    """

    def __init__(self, cal_obj, parent_job_dir='.',
                 job_dir='./MeasurementSolvation',
                 sol_params=None):
        Measurement.__init__(self, cal_objs=cal_obj,
                             parent_job_dir=parent_job_dir,
                             job_dir=job_dir)
        self.sol_params = sol_params or {'EB_K': [78.4],
                                         'TAU': [0],
                                         'LAMBDA_D_K': [3.0],
                                         'NELECT': []}

    def setup(self):
        """
        setup solvation jobs for the calibrate objects
        copies WAVECAR and sets the solvation params in the incar file
        also dumps system.json file in each directory for the database
        crawler
        mind: works only for cal objects that does only single
        calculations
        """
        for cal in self.cal_objs:
            jdir = cal.old_job_dir_list[0]
            cal.poscar = Poscar.from_file(jdir + os.sep + 'POSCAR')
            cal.potcar = Potcar.from_file(jdir + os.sep + 'POTCAR')
            cal.kpoints = Kpoints.from_file(jdir + os.sep + 'KPOINTS')
            cal.incar = Incar.from_file(jdir + os.sep + 'INCAR')
            cal.incar['LSOL'] = '.TRUE.'
            syms = [site.specie.symbol for site in cal.poscar.structure]
            zvals = {p.symbol: p.nelectrons for p in cal.potcar}
            nelectrons = sum([zvals[a[0]] * len(tuple(a[1]))
                              for a in itertools.groupby(syms)])
            keys = [k for k in self.sol_params.keys()
                    if self.sol_params[k]]
            prod_list = [self.sol_params.get(k) for k in keys]
            for params in itertools.product(*tuple(prod_list)):
                job_dir = self.job_dir + os.sep \
                    + cal.old_job_dir_list[0].replace(os.sep,
                                                      '_').replace('.',
                                                                   '_') \
                    + os.sep + 'SOL'
                for i, k in enumerate(keys):
                    if k == 'NELECT':
                        cal.incar[k] = params[i] + nelectrons
                    else:
                        cal.incar[k] = params[i]
                    job_dir = job_dir + os.sep + k + os.sep + str(
                        cal.incar[k]).replace('.', '_')
                if not os.path.exists(job_dir):
                    os.makedirs(job_dir)
                with open(job_dir + os.sep + 'system.json', 'w') as f:
                    json.dump(dict(list(zip(keys, params))), f)
                wavecar_file = cal.old_job_dir_list[0] + os.sep + 'WAVECAR'
                if os.path.isfile(wavecar_file):
                    shutil.copy(wavecar_file, job_dir + os.sep + 'WAVECAR')
                    cal.add_job(job_dir=job_dir)
                else:
                    logger.critical('WAVECAR doesnt exist. Aborting ...')
                    sys.exit(0)

    def make_measurements(self):
        """
        get solvation energies
        """
        pass


class MeasurementInterface(Measurement):
    """
    Interface

    Takes list of Calibration Objects of Interface, Slab and
    Ligand and separates them

    Args:
        cal_objs: List of Calibration Objects

    """

    def __init__(self, cal_objs, parent_job_dir='.',
                 job_dir='./MeasurementInterface'):
        Measurement.__init__(self, cal_objs=cal_objs,
                             parent_job_dir=parent_job_dir,
                             job_dir=job_dir)
        self.cal_slabs = []
        self.cal_interfaces = []
        self.cal_ligands = []
        for cal in self.cal_objs:
            if isinstance(cal, CalibrateSlab):
                self.cal_slabs.append(cal)
            elif isinstance(cal, CalibrateInterface):
                self.cal_interfaces.append(cal)
            elif isinstance(cal, CalibrateMolecule):
                self.cal_ligands.append(cal)

    def setup(self):
        """
        setup static jobs for the calibrate objects
        copies CONTCAR to POSCAR
        sets NSW = 0
        write system.json file for database crawler
        """
        d = {}
        for cal in self.cal_objs:
            for i, jdir in enumerate(cal.old_job_dir_list):
                job_dir = self.job_dir + os.sep \
                    + jdir.replace(os.sep, '_').replace('.', '_') + \
                    os.sep + 'STATIC'
                cal.incar = Incar.from_file(jdir + os.sep + 'INCAR')
                cal.incar['EDIFF'] = '1E-6'
                cal.incar['NSW'] = 0
                cal.potcar = Potcar.from_file(jdir + os.sep + 'POTCAR')
                cal.kpoints = Kpoints.from_file(jdir + os.sep + 'KPOINTS')
                contcar_file = jdir + os.sep + 'CONTCAR'
                if os.path.isfile(contcar_file):
                    cal.poscar = Poscar.from_file(contcar_file)
                    if cal in self.cal_slabs or cal in self.cal_interfaces:
                        try:
                            d['hkl'] = cal.system['hkl']
                        except:
                            logger.critical("""the calibrate object
                            doesnt have a system set for calibrating""")
                    if cal in self.cal_interfaces:
                        try:
                            d['ligand'] = cal.system['ligand']['name']
                        except:
                            logger.critical("""the calibrate object
                            doesnt have a system set for calibrating""")
                    if not os.path.exists(job_dir):
                        os.makedirs(job_dir)
                    if d:
                        with open(job_dir + os.sep + 'system.json', 'w') as f:
                            json.dump(d, f)
                    cal.add_job(job_dir=job_dir)
                else:
                    logger.critical("""CONTCAR doesnt exist.
                    Setting up job using input set in the old
                    calibration directory""")
                    cal.poscar = Poscar.from_file(jdir + os.sep + 'POSCAR')
                    cal.add_job(job_dir=job_dir)

    def make_measurements(self):
        """
        get the slab, ligand and interface energies
        compute binding energies
        """
        E_interfaces = {}
        E_slabs = {}
        E_ligands = {}
        E_binding = {}
        for cal in self.cal_slabs:
            key = str(cal.system['hkl'])
            E_slabs[key] = self.get_energy(cal)
        for cal in self.cal_ligands:
            key = cal.system['ligand']['name']
            E_ligands[key] = self.get_energy(cal)
        for cal in self.cal_interfaces:
            key_slab = str(cal.system['hkl'])
            key_ligand = cal.system['ligand']['name']
            key = key_slab + key_ligand
            E_interfaces[key] = self.get_energy(cal)
            E_binding[key] = E_interfaces[key] \
                - E_slabs[key_slab] \
                - cal.system['num_ligands'] * E_ligands[
                key_ligand]
        logger.info('Binding energy = {}'.format(E_binding))


# test
if __name__ == '__main__':
    from pymatgen.core.structure import Structure, Molecule
    from mpinterfaces.interface import Ligand

    # PbS 100 surface with single hydrazine as ligand
    strt = Structure.from_file("POSCAR.mp-21276_PbS")
    mol_struct = Structure.from_file("POSCAR_diacetate")
    mol = Molecule(mol_struct.species, mol_struct.cart_coords)
    hydrazine = Ligand([mol])
    supercell = [1, 1, 1]
    hkl = [1, 1, 1]
    min_thick = 10
    min_vac = 12
    surface_coverage = 0.01
    adsorb_on_species = 'S'
    adatom_on_lig = 'Pb'
    displacement = 3.0
    iface = Interface(strt, hkl=hkl, min_thick=min_thick,
                      min_vac=min_vac,
                      supercell=supercell, surface_coverage=0.01,
                      ligand=hydrazine, displacement=displacement,
                      adatom_on_lig=adatom_on_lig,
                      adsorb_on_species=adsorb_on_species,
                      primitive=False, coverage_tol=0.5)
    iface.create_interface()
    iface.sort()

    incarparams = {'System': 'test',
                   'ENCUT': 400,
                   'ISMEAR': 1,
                   'SIGMA': 0.1,
                   'EDIFF': 1E-6}
    incar = Incar(params=incarparams)
    poscar = Poscar(iface)
    potcar = Potcar(symbols=poscar.site_symbols, functional='PBE',
                    sym_potcar_map=None)
    kpoints = Kpoints.monkhorst_automatic(kpts=(16, 16, 16),
                                          shift=(0, 0, 0))
    cal = CalibrateInterface(incar, poscar, potcar, kpoints,
                             system=iface.as_dict(),
                             job_dir='test', job_cmd=['ls', '-lt'])
    cal.setup()
    cal.run()
    # list of calibrate objects
    cal_objs = [cal]
    # check whether the cal jobs were done
    # Calibrate.check_calcs(cal_objs)
    # set the measurement
    # measure = MeasurementInterface(cal_objs, job_dir='./Measurements')
    measure = MeasurementSolvation(cal_objs, job_dir='./MSR_SOL',
                                   sol_params={'EB_K': [78.4],
                                               'TAU': [0],
                                               'LAMBDA_D_K': [3.0],
                                               'NELECT': [1, -1]})
    measure.setup()
    measure.run()
