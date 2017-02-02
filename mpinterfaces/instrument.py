# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
The instrument module:
defines the inputset and the job

"""

import sys
import os, shutil
import subprocess
import logging

from pymatgen.io.vasp.inputs import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import DictSet

from custodian.custodian import Job, ErrorHandler

from monty.json import MontyDecoder

from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from mpinterfaces.data_processor import MPINTVasprun

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


class MPINTVaspInputSet(DictSet):
    """
    defines the set of input required for a vasp job i.e
    create INCAR, POSCAR, POTCAR & KPOINTS files
    """

    def __init__(self, name, incar, poscar, potcar, kpoints,
                 qadapter=None, script_name='submit_script',
                 vis_logger=None, reuse_path=None, **kwargs):
        """
        default INCAR from config_dict

        """
        self.name = name
        self.incar = Incar.from_dict(incar.as_dict())
        self.poscar = Poscar.from_dict(poscar.as_dict())
        self.potcar = Potcar.from_dict(potcar.as_dict())
        if not type(kpoints) == str:
           self.kpoints = Kpoints.from_dict(kpoints.as_dict())
        else:
           self.kpoints = kpoints
        self.reuse_path = reuse_path # complete reuse paths
        self.extra = kwargs
        if qadapter is not None:
            self.qadapter = qadapter.from_dict(qadapter.to_dict())
        else:
            self.qadapter = None
        self.script_name = script_name
        config_dict = {}
        config_dict['INCAR'] = self.incar.as_dict()
        config_dict['POSCAR'] = self.poscar.as_dict()
        # caution the key and the value are not always the same
        config_dict['POTCAR'] = self.potcar.as_dict()
        # dict(zip(self.potcar.as_dict()['symbols'],
        # self.potcar.as_dict()['symbols']))
        if not type(kpoints) == str:
            config_dict['KPOINTS'] = self.kpoints.as_dict()
        else:
            config_dict['KPOINTS'] = self.kpoints
        # self.user_incar_settings = self.incar.as_dict()
        DictSet.__init__(self, name, config_dict,
                                  ediff_per_atom=False, **kwargs)
        if vis_logger:
            self.logger = vis_logger
        else:
            self.logger = logger

    def write_input(self, job_dir, make_dir_if_not_present=True,
                    write_cif=False):
        """
        the input files are written to the job_dir
        process(if needed) and write the input files in each directory
        structures read from the poscar files in the directory
        """
        d = job_dir
        if make_dir_if_not_present and not os.path.exists(d):
            os.makedirs(d)
        self.logger.info('writing inputset to : ' + d)
        self.incar.write_file(os.path.join(d, 'INCAR'))

        if not type(self.kpoints) == str:
            ## maybe temporary fix, pymatgen does not seem
            ## to have a versatile kpoints object for writing a
            ## HSE Kpoints file
            self.kpoints.write_file(os.path.join(d, 'KPOINTS'))
        else:
            with open(os.path.join(d, 'KPOINTS'), 'w') as kpts:
                for line in self.kpoints:
                   kpts.write(line)

        self.potcar.write_file(os.path.join(d, 'POTCAR'))
        self.poscar.write_file(os.path.join(d, 'POSCAR'),
                               significant_figures=10)
        if self.qadapter is not None:
            with open(os.path.join(d, self.script_name), 'w') as f:
                queue_script = self.qadapter.get_script_str(job_dir)
                f.write(queue_script)

    def as_dict(self):
        qadapter = None
        if self.qadapter:
            qadapter = self.qadapter.to_dict()

        if not type(self.kpoints) == str:
           kpoints = self.kpoints.as_dict()
        else:
           kpoints = [self.kpoints]

        d = dict(name=self.name, incar=self.incar.as_dict(),
                 poscar=self.poscar.as_dict(),
                 potcar=self.potcar.as_dict(),
                 kpoints=kpoints,
                 qadapter=qadapter, script_name=self.script_name,
                 kwargs=self.extra)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["logger"] = self.logger.name
        return d

    @classmethod
    def from_dict(cls, d):
        incar = Incar.from_dict(d["incar"])
        poscar = Poscar.from_dict(d["poscar"])
        potcar = Potcar.from_dict(d["potcar"])
        kpoints = Kpoints.from_dict(d["kpoints"])
        qadapter = None
        if d["qadapter"] is not None:
            qadapter = CommonAdapter.from_dict(d["qadapter"])
        script_name = d["script_name"]
        return MPINTVaspInputSet(d["name"], incar, poscar, potcar,
                                 kpoints, qadapter,
                                 script_name=script_name,
                                 vis_logger=logging.getLogger(d["logger"]),
                                 **d["kwargs"])


class MPINTJob(Job):
    """
    defines a job i.e setup the required input files and
    launch the job

    Args:
       job_cmd: a list, the command to be issued in each job_dir
                 eg: ['qsub', 'submit_job']
       job_dir: the directory from which the jobs will be launched
    """

    def __init__(self, job_cmd, name='noname', output_file="job.out",
                 parent_job_dir='.', job_dir='untitled', suffix="",
                 final=True, gzipped=False, backup=False, vis=None,
                 auto_npar=True, settings_override=None, wait=True,
                 vjob_logger=None):
        self.job_cmd = job_cmd
        self.name = name
        self.output_file = output_file
        self.parent_job_dir = parent_job_dir
        self.job_dir = job_dir
        self.final = final
        self.backup = backup
        self.gzipped = gzipped
        self.vis = vis
        self.suffix = suffix
        self.settings_override = settings_override
        self.auto_npar = auto_npar
        self.wait = wait
        if vjob_logger:
            self.logger = vjob_logger
        else:
            self.logger = logger

    def setup(self):
        """
        write the input files to the job_dir
        """
        self.vis.write_input(self.job_dir)
        if self.backup:
            os.chdir(os.path.abspath(self.job_dir))
            for f in os.listdir('.'):
                shutil.copy(f, "{}.orig".format(f))
            os.chdir(self.parent_job_dir)

    def run(self):
        """
        move to the job_dir, launch the job and back to the
        parent job directory
        """
        os.chdir(os.path.abspath(self.job_dir))
        self.logger.info('running in : ' + self.job_dir)
        p = None
        # if launching jobs via batch system
        if self.vis.qadapter is not None:
            submit_cmd = \
                self.vis.qadapter.q_commands[self.vis.qadapter.q_type][
                    "submit_cmd"]
            cmd = [submit_cmd, self.vis.script_name]
            with open(self.output_file, 'w') as f:
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
                stdout, stderr = p.communicate()
                self.job_id = stdout.rstrip('\n').split()[-1]
                f.write(self.job_id)
        else:
            cmd = list(self.job_cmd)
            with open(self.output_file, 'w') as f:
                p = subprocess.Popen(cmd, stdout=f, stderr=f)
            self.job_id = 0  # None
        os.chdir(self.parent_job_dir)
        if self.wait:
            return p
        else:
            return 0

    def postprocess(self):
        pass

    def name(self):
        return self.__class__.__name__

    def as_dict(self):
        d = dict(job_cmd=self.job_cmd, name=self.name,
                 output_file=self.output_file,
                 parent_job_dir=self.parent_job_dir,
                 job_dir=self.job_dir, suffix=self.suffix,
                 final=self.final, gzipped=self.gzipped,
                 backup=self.backup, vis=self.vis.as_dict(),
                 auto_npar=self.auto_npar,
                 settings_override=self.settings_override,
                 wait=self.wait)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["logger"] = self.logger.name
        return d

    @classmethod
    def from_dict(cls, d):
        vis = MontyDecoder().process_decoded(d["vis"])
        return MPINTVaspJob(d["job_cmd"], name=d["name"],
                            output_file=d["output_file"],
                            parent_job_dir=d["parent_job_dir"],
                            job_dir=d["job_dir"],
                            suffix=d["suffix"], final=d["final"],
                            gzipped=d["gzipped"],
                            backup=d["backup"], vis=vis,
                            auto_npar=d["auto_npar"],
                            settings_override=d["settings_override"],
                            wait=d["wait"],
                            vjob_logger=logging.getLogger(d["logger"]))


class MPINTVaspJob(MPINTJob):
    """
    defines a vasp job i.e setup the required input files and
    launch the job

    Args:
       job_cmd: a list, the command to be issued in each job_dir
                 eg: ['qsub', 'submit_job']
       job_dir: the directory from which the jobs will be launched
    """

    def __init__(self, job_cmd, name='noname', output_file="job.out",
                 parent_job_dir='.', job_dir='untitled', suffix="",
                 final=True, gzipped=False, backup=False, vis=None,
                 auto_npar=True, settings_override=None, wait=True,
                 vjob_logger=None):
        MPINTJob.__init__(self, job_cmd, name=name,
                          output_file=output_file,
                          parent_job_dir=parent_job_dir,
                          job_dir=job_dir, suffix=suffix,
                          final=final, gzipped=gzipped,
                          backup=backup, vis=vis, auto_npar=auto_npar,
                          settings_override=settings_override,
                          wait=wait, vjob_logger=vjob_logger)

    def get_final_energy(self):
        vasprun_file_path = self.job_dir + os.sep + 'vasprun.xml'
        try:
            vasprun = MPINTVasprun(vasprun_file_path,
                                   parse_potcar_file=False)
            if vasprun.converged:
                self.logger.info("job {0} in {1} converged".format(self.job_id,
                                                                   self.job_dir))
                return vasprun.final_energy
            else:
                self.logger.info(
                    "job {0} in {1} NOT converged".format(self.job_id,
                                                          self.job_dir))
                return None
        except Exception as ex:
            self.logger.info(
                "error reading vasprun.xml, probably the job {0} in {1} is not done yet.".format(
                    self.job_id,
                    self.job_dir))
            return None


class MPINTVaspErrors(ErrorHandler):
    """
    handles restarting of jobs that exceed the walltime
    employs the check + correct method of custodian ErrorHandler
    """
    pass
