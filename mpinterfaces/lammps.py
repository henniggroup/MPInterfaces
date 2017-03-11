# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Calibrate LAMMPS jobs
"""

from six.moves import map
from six.moves import zip

import os
import logging
from collections import OrderedDict

from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from monty.json import MSONable, MontyDecoder

from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from ase.calculators.lammpsrun import LAMMPS, prism

from mpinterfaces.instrument import MPINTJob
from mpinterfaces.calibrate import Calibrate

from mpinterfaces.default_logger import get_default_logger

logger = get_default_logger(__name__)


class MPINTLammps(LAMMPS, MSONable):
    """
    setup LAMMPS for given structure and parameters
    extends ase.calculators.lammpsrun.LAMMPS
    """

    def __init__(self, structure, parameters={},
                 label='mpintlmp', specorder=None,
                 always_triclinic=False, no_data_file=False):
        LAMMPS.__init__(self, label=label,
                        parameters=parameters,
                        specorder=specorder, files=[],
                        always_triclinic=always_triclinic,
                        no_data_file=no_data_file)
        self.structure = structure
        self.atoms = AseAtomsAdaptor().get_atoms(structure)
        self.label = label
        self.parameters = parameters
        self.specorder = specorder
        self.always_triclinic = always_triclinic
        self.no_data_file = no_data_file
        self.charges = None
        if 'charges' in self.parameters:
            self.charges = self.parameters['charges']

    def write_lammps_data(self, f):
        """
        write lammps structure data
        from ase with custom modifications
        """
        f.write(f.name + ' (written by mpinterfaces) \n\n')
        symbols = self.atoms.get_chemical_symbols()
        n_atoms = len(symbols)
        f.write('%d \t atoms \n' % n_atoms)
        if self.specorder is None:
            species = [tos.symbol
                       for tos in self.structure.types_of_specie]
        else:
            species = self.specorder
        n_atom_types = len(species)
        f.write('%d  atom types\n' % n_atom_types)
        p = prism(self.atoms.get_cell())
        xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism_str()
        f.write('0.0 %s  xlo xhi\n' % xhi)
        f.write('0.0 %s  ylo yhi\n' % yhi)
        f.write('0.0 %s  zlo zhi\n' % zhi)
        if self.always_triclinic or p.is_skewed():
            f.write('%s %s %s  xy xz yz\n' % (xy, xz, yz))
        f.write('\n\n')
        f.write('Atoms \n\n')
        for i, r in enumerate(map(p.pos_to_lammps_str,
                                  self.atoms.get_positions())):
            c = 0.0
            if self.charges:
                c = self.charges[symbols[i]]
            s = species.index(symbols[i]) + 1
            if 'atom_style' in self.parameters:
                if self.parameters['atom_style'] == 'charge':

                    f.write('%6d %3d %6f %s %s %s\n' %
                            ((i + 1, s, c) + tuple(r)))
                else:
                    f.write('%6d %3d %s %s %s\n' % ((i + 1, s) + tuple(r)))
        if self.atoms.get_velocities() is not None:
            f.write('\n\nVelocities \n\n')
            for i, v in enumerate(self.atoms.get_velocities()):
                f.write('%6d %s %s %s\n' % ((i + 1,) + tuple(v)))
        f.close()

    def write_lammps_in(self, lammps_in=None, lammps_trj=None,
                        lammps_data=None):
        """
        write lammps input file
        from ase with custom modifications
        """
        f = lammps_in
        f.write('clear\n' +
                ('variable dump_file string "%s"\n' % lammps_trj) +
                ('variable data_file string "%s"\n' % lammps_data))
        parameters = self.parameters
        pbc = self.atoms.get_pbc()
        f.write('units metal \n')
        if 'atom_style' in parameters:
            f.write('atom_style %s \n' % parameters['atom_style'])
        else:
            f.write('atom_style atomic \n')
        if 'boundary' in parameters:
            f.write('boundary %s \n' % parameters['boundary'])
        else:
            f.write('boundary %c %c %c \n' %
                    tuple('sp'[x] for x in pbc))
        f.write('atom_modify sort 0 0.0 \n')
        for key in ('neighbor', 'newton'):
            if key in parameters:
                f.write('%s %s \n' % (key, parameters[key]))
        f.write('\n')
        # If no_data_file,
        # write the simulation box and the atoms
        if self.no_data_file:
            p = self.prism
            f.write('lattice sc 1.0\n')
            xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism_str()
            if self.always_triclinic or p.is_skewed():
                f.write('region asecell prism 0.0 %s 0.0 %s 0.0 %s ' %
                        (xhi, yhi, zhi))
                f.write('%s %s %s side in units box\n' % (xy, xz, yz))
            else:
                f.write(('region asecell block 0.0 %s 0.0 %s 0.0 %s '
                         'side in units box\n') % (xhi, yhi, zhi))
            symbols = self.atoms.get_chemical_symbols()
            if self.specorder is None:
                species = [tos.symbol for tos in
                           self.structure.types_of_specie]
            else:
                species = self.specorder
            n_atom_types = len(species)
            species_i = dict([(s, i + 1) for i, s in enumerate(species)])
            f.write('create_box %i asecell\n' % n_atom_types)
            for s, pos in zip(symbols, self.atoms.get_positions()):
                f.write('create_atoms %i single %s %s %s units box\n' %
                        ((species_i[s],) + p.pos_to_lammps_fold_str(pos)))
        else:
            f.write('read_data %s\n' % lammps_data)
        # interaction
        f.write('\n### interactions \n')
        if (('pair_style' in parameters) and ('pair_coeff' in parameters)):
            pair_style = parameters['pair_style']
            f.write('pair_style %s \n' % pair_style)
            for pair_coeff in parameters['pair_coeff']:
                f.write('pair_coeff %s \n' % pair_coeff)
            if 'mass' in parameters:
                for mass in parameters['mass']:
                    f.write('mass %s \n' % mass)
        else:
            # default interaction
            f.write('pair_style lj/cut 2.5 \n' +
                    'pair_coeff * * 1 1 \n' +
                    'mass * 1.0 \n')
        f.write('\n### run\n')
        if 'fix' in parameters:
            if parameters['fix']:
                for i in parameters['fix']:
                    f.write('fix %s\n' % i)
        else:
            f.write('fix fix_nve all nve\n')
        if 'thermo_style' in parameters:
            f.write('thermo_style %s\n' % parameters['thermo_style'])
        else:
            f.write(('thermo_style custom %s\n') %
                    (' '.join(self._custom_thermo_args)))
        if 'thermo_modify' in parameters:
            f.write('thermo_modify %s\n' % parameters['thermo_modify'])
        else:
            f.write('thermo_modify flush yes\n')
        if 'thermo' in parameters:
            f.write('thermo %s\n' % parameters['thermo'])
        else:
            f.write('thermo 1\n')
        if 'minimize' in parameters:
            f.write('minimize %s\n' % parameters['minimize'])
        if 'run' in parameters:
            f.write('run %s\n' % parameters['run'])
        if not (('minimize' in parameters) or ('run' in parameters)):
            f.write('run 0\n')
        if 'dump' in parameters:
            f.write('dump %s\n' % parameters['dump'])
        else:
            f.write(
                'dump dump_all all custom 1 %s id type x y z vx vy vz fx fy fz\n' % lammps_trj)
        f.write('print __end_of_ase_invoked_calculation__\n')
        f.write('log /dev/stdout\n')
        f.close()

    def as_dict(self):
        d = dict(structure=self.structure.as_dict(),
                 parameters=self.parameters, label=self.label,
                 specorder=self.specorder,
                 always_triclinic=self.always_triclinic,
                 no_data_file=self.no_data_file)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        structure = Structure.from_dict(d["structure"])
        return MPINTLammps(structure, parameters=d["parameters"],
                           label=d["label"],
                           specorder=d["specorder"],
                           always_triclinic=d["always_triclinic"],
                           no_data_file=d["no_data_file"])


class MPINTLammpsInput(MSONable):
    """
    create inputs(data and input file) for LAMMPS
    """

    def __init__(self, mplmp, qadapter=None, vis_logger=None):
        self.mplmp = mplmp
        if qadapter is not None:
            self.qadapter = qadapter.from_dict(qadapter.to_dict())
        else:
            self.qadapter = None
        if vis_logger:
            self.logger = vis_logger
        else:
            self.logger = logger

    def write_input(self, job_dir, make_dir_if_not_present=True):
        """
        write LAMMPS input set to job_dir
        """
        d = job_dir
        if make_dir_if_not_present and not os.path.exists(d):
            os.makedirs(d)
        self.logger.info('writing inputset to : ' + d)
        lmp_in = os.path.join(d, 'inp')
        lmp_data_name = 'data'
        lmp_data = os.path.join(d, lmp_data_name)
        lmp_trj_name = 'trj'
        lmp_in_fd = open(lmp_in, 'w')
        lmp_data_fd = open(lmp_data, 'w')
        self.mplmp.write_lammps_data(lmp_data_fd)
        self.mplmp.write_lammps_in(lammps_in=lmp_in_fd,
                                   lammps_trj=lmp_trj_name,
                                   lammps_data=lmp_data_name)
        if self.qadapter is not None:
            self.script_name = 'submit_script'
            with open(os.path.join(d, self.script_name), 'w') as f:
                queue_script = self.qadapter.get_script_str(job_dir)
                f.write(queue_script)
        lmp_in_fd.close()
        lmp_data_fd.close()

    def as_dict(self):
        qadapter = None
        if self.qadapter:
            qadapter = self.qadapter.to_dict()
        d = dict(qadapter=qadapter, mplmp=self.mplmp.as_dict())
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["logger"] = self.logger.name
        return d

    @classmethod
    def from_dict(cls, d):
        qadapter = None
        if d["qadapter"] is not None:
            qadapter = CommonAdapter.from_dict(d["qadapter"])
        return MPINTLammpsInput(
            MPINTLammps.from_dict(d["mplmp"]),
            qadapter,
            vis_logger=logging.getLogger(d["logger"]))


class MPINTLammpsJob(MPINTJob):
    """
    Define LAMMPS job
    """

    def __init__(self, job_cmd, name='mpintlmpjob',
                 output_file="job.out", parent_job_dir='.',
                 job_dir='untitled', final=True, gzipped=False,
                 backup=False, vis=None, settings_override=None,
                 wait=True, vjob_logger=None):
        MPINTJob.__init__(self, job_cmd, name=name,
                          output_file=output_file,
                          parent_job_dir=parent_job_dir,
                          job_dir=job_dir,
                          final=final, gzipped=gzipped,
                          backup=backup, vis=vis,
                          settings_override=settings_override,
                          wait=wait, vjob_logger=vjob_logger)

    def get_final_energy(self, lammps_log='log.lammps'):
        """
        return the final total energy
        """
        energy = None
        try:
            f = open(os.path.join(self.job_dir, lammps_log), 'r')
            self.vis.mplmp.read_lammps_log(lammps_log=f,
                                           PotEng_first=False)
            energy = self.vis.mplmp.thermo_content[-1]['etotal']
            f.close()
        except:
            energy = None
        return energy

    def as_dict(self):
        d = dict(job_cmd=self.job_cmd,
                 output_file=self.output_file,
                 parent_job_dir=self.parent_job_dir,
                 job_dir=self.job_dir, final=self.final,
                 gzipped=self.gzipped, backup=self.backup,
                 vis=self.vis.as_dict(),
                 settings_override=self.settings_override,
                 wait=self.wait)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["logger"] = self.logger.name
        return d

    @classmethod
    def from_dict(cls, d):
        vis = MontyDecoder().process_decoded(d["vis"])
        return MPINTLammpsJob(d["job_cmd"],
                              output_file=d["output_file"],
                              parent_job_dir=d["parent_job_dir"],
                              job_dir=d["job_dir"],
                              final=d["final"],
                              gzipped=d["gzipped"],
                              backup=d["backup"], vis=vis,
                              settings_override=d["settings_override"],
                              wait=d["wait"],
                              vjob_logger=logging.getLogger(
                                  d["logger"]))


class CalibrateLammps(Calibrate):
    """
    Defines LAMMPS workflow consisting of LAMMPS jobs
    """

    def __init__(self, parameters, structure=None, parent_job_dir='.',
                 job_dir='./cal_lammps', qadapter=None,
                 job_cmd='qsub', wait=True, is_matrix=False,
                 turn_knobs=OrderedDict([('STRUCTURES', []),
                                         ('PAIR_COEFF', [])]),
                 checkpoint_file=None, cal_logger=None):
        self.parameters = parameters
        self.structure = structure
        Calibrate.__init__(self, None, None, None, None,
                           parent_job_dir=parent_job_dir,
                           job_dir=job_dir, qadapter=qadapter,
                           job_cmd=job_cmd, wait=wait,
                           turn_knobs=turn_knobs,
                           is_matrix=is_matrix,
                           checkpoint_file=checkpoint_file,
                           cal_logger=cal_logger)

    def val_to_name(self, val):
        if isinstance(val, Structure):
            return '_'.join([val.formula.replace(' ', ''),
                             str(val.volume)])
        else:
            name = str(val).replace(' ', '_')
            return name.split('/')[-1]

    def key_to_name(self, key):
        return str(key)

    def setup_structure_jobs(self, structures, paircoeff):
        for s in structures:
            self.structure = s
            self.set_paircoeff(s, paircoeff)
            if not self.is_matrix:
                name = self.val_to_name(s)
                job_dir = self.job_dir + os.sep + \
                    self.key_to_name('STRUCTURES') + os.sep + name
                self.add_job(name=name, job_dir=job_dir)

    def setup_paircoeff_jobs(self, paircoeff_list):
        for pcoeff in paircoeff_list:
            self.set_paircoeff(self.structure, pcoeff)
            if not self.is_matrix:
                name = self.val_to_name(pcoeff)
                job_dir = self.job_dir + os.sep + \
                    self.key_to_name('PAIR_COEFF') + os.sep + name
                self.add_job(name=name, job_dir=job_dir)

    def setup_params_jobs(self, params_list):
        for params in params_list:
            self.parameters = params
            if not self.is_matrix:
                name = self.val_to_name(params)
                job_dir = self.job_dir + os.sep + \
                    self.key_to_name('PARAMS') + os.sep + name
                self.add_job(name=name, job_dir=job_dir)

    def setup_genericparam_jobs(self, key, vals):
        for v in vals:
            self.parameters[key.lower()] = v
            if not self.is_matrix:
                name = self.val_to_name(v)
                job_dir = self.job_dir + os.sep + \
                    self.key_to_name(key) + os.sep + name
                self.add_job(name=name, job_dir=job_dir)

    def set_paircoeff(self, structure, pcoeff):
        types_of_species = ' '.join(
            [tos.symbol for tos in structure.types_of_specie])
        atomic_mass = [
            str(i + 1) + ' ' + tos.atomic_mass.__repr__()
            for i, tos in enumerate(structure.types_of_specie)]
        self.parameters['mass'] = atomic_mass
        self.parameters['pair_coeff'] = [
            '* * {0} {1}'.format(pcoeff, types_of_species)]

    def add_job(self, name='noname', job_dir='.'):
        """
        add lammps job given MPINTLammps object
        """
        mplmp = MPINTLammps(self.structure,
                            parameters=self.parameters)
        lmp_inp = MPINTLammpsInput(mplmp, self.qadapter, self.logger)
        job = MPINTLammpsJob(self.job_cmd, name=name, final=True,
                             parent_job_dir=self.parent_job_dir,
                             job_dir=job_dir, vis=lmp_inp,
                             wait=self.wait, vjob_logger=self.logger)
        self.jobs.append(job)

    def _setup(self, turn_knobs=None):
        """
        setup workflow for the given turn_knobs i.e create jobs
        for each knob settings
        keys with special support(everything else handled by
        genericparam job setup):
             STRUCTURES: list of pymatgen structure objects
             PARAMS: list of lammps parameter dictionaries
             PAIR_COEFF: list of pair coefficient files
        """
        if turn_knobs is None:
            turn_knobs = self.turn_knobs
        if any(turn_knobs.values()):
            for k, v in turn_knobs.items():
                if k == 'STRUCTURES' and v:
                    self.setup_structure_jobs(
                        v, self.turn_knobs['PAIR_COEFF'][0])
                elif k == 'PAIR_COEFF' and v:
                    self.setup_paircoeff_jobs(v)
                elif k == 'PARAMS' and v:
                    self.setup_params_jobs(v)
                else:
                    self.setup_genericparam_jobs(k, v)
        else:
            self.logger.warn('knobs not set, running a single job')
            self.add_job(name='single_job', job_dir=self.job_dir)

    def as_dict(self):
        qadapter = None
        system = None
        if self.qadapter:
            qadapter = self.qadapter.to_dict()
        if self.system is not None:
            system = self.system
        d = dict(
            parent_job_dir=self.parent_job_dir,
            job_dir=self.job_dir,
            qadapter=qadapter, job_cmd=self.job_cmd,
            wait=self.wait,
            turn_knobs=self.turn_knobs,
            job_ids=self.job_ids)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        cal = CalibrateLammps(parent_job_dir=d["parent_job_dir"],
                              job_dir=d["job_dir"],
                              qadapter=d.get("qadapter"),
                              job_cmd=d["job_cmd"], wait=d["wait"],
                              turn_knobs=d["turn_knobs"])
        cal.job_ids = d["job_ids"]
        return cal
