# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Defines various firetasks
"""

from pymatgen.apps.borg.queen import BorgQueen

from monty.json import MontyDecoder

from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import FWSerializable
from fireworks.utilities.fw_utilities import explicit_serialize

from mpinterfaces.database import MPINTVaspToDbTaskDrone
from mpinterfaces.default_logger import get_default_logger

logger = get_default_logger(__name__)


def load_class(mod, name):
    """
    load class named name from module mod
    this function is adapted from materialsproject
    """
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)


def get_cal_obj(d):
    """
    construct a calibration object from the input dictionary, d

    returns a calibration object
    """
    cal = MontyDecoder().process_decoded(d)
    # default
    # if not d.get("qadapter"):
    #    qadapter, job_cmd = get_run_cmmnd(**d.get('que_params', {}))
    #    cal.qadapter = qadapter
    #    cal.job_cmd = job_cmd
    return cal


@explicit_serialize
class MPINTCalibrateTask(FireTaskBase):
    """
    Calibration Task
    """

    optional_params = ["que_params"]

    def run_task(self, fw_spec):
        """
        launch jobs to the queue
        """
        cal = get_cal_obj(self)
        cal.setup()
        cal.run()
        d = cal.as_dict()
        d.update({'que_params': self.get('que_params')})
        return FWAction(mod_spec=[{'_push': {'cal_objs': d}}])


@explicit_serialize
class MPINTMeasurementTask(FireTaskBase, FWSerializable):
    """
    Measurement Task
    """
    required_params = ["measurement"]
    optional_params = ["que_params", "job_cmd", "other_params", "fw_id"]

    def run_task(self, fw_spec):
        """
        setup up a measurement task using the prior calibration jobs
        and run
        """
        cal_objs = []
        logger.info(
            'The measurement task will be constructed from {} calibration objects'
            .format(len(fw_spec['cal_objs'])))
        for calparams in fw_spec['cal_objs']:
            calparams.update({'que_params': self.get('que_params')})
            cal = get_cal_obj(calparams)
            cal_objs.append(cal)
        done = load_class("mpinterfaces.calibrate", "Calibrate").check_calcs(
            cal_objs)
        if not done:
            logger.info('Calibration not done yet. Try again later')
            logger.info('All subsequent fireworks will be defused')
            logger.info("""Try re-running this firework again later.
            Re-running this firework will activate all the subsequent foreworks too""")
            logger.info('This fireworks id = {}'.format(self.get("fw_id")))
            return FWAction(defuse_children=True)
            # to enable dynamic workflow, uncomment the following
            # if self.get("fw_id"):
            #    fw_id = int(self.get("fw_id")) + 1
            #    self["fw_id"] = fw_id
            #    new_fw = Firework(MPINTMeasurementTask(self),
            #                      spec={'cal_objs':fw_spec['cal_objs']},
            #                      name = 'new_fw', fw_id = -fw_id)
            # else:
            #    new_fw = Firework(MPINTMeasurementTask(self),
            #                      spec={'cal_objs':fw_spec['cal_objs']},
            #                      name = 'new_fw')
            #
            # return FWAction(detours=new_fw)
        else:
            measure = load_class("mpinterfaces.measurement",
                                 self['measurement'])(cal_objs, **self.get(
                                     "other_params", {}))
            job_cmd = None
            if self.get("job_cmd", None) is not None:
                job_cmd = self.get("job_cmd")
            measure.setup()
            measure.run(job_cmd=job_cmd)
            cal_list = []
            for cal in measure.cal_objs:
                d = cal.as_dict()
                d.update({'que_params': self.get('que_params')})
                cal_list.append(d)
            return FWAction(update_spec={'cal_objs': cal_list})


@explicit_serialize
class MPINTDatabaseTask(FireTaskBase, FWSerializable):
    """
    submit data to the database firetask
    """
    required_params = ["measure_dir"]
    optional_params = ["dbase_params"]

    def run_task(self, fw_spec):
        """
        go through the measurement job dirs and
        put the measurement jobs in the database
        """
        drone = MPINTVaspToDbTaskDrone(**self.get("dbase_params", {}))
        queen = BorgQueen(drone)  # , number_of_drones=ncpus)
        queen.serial_assimilate(self["measure_dir"])
        return FWAction()
