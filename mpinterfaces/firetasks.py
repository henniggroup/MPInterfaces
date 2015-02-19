from __future__ import division, unicode_literals, print_function

"""
Defines various firetasks
"""

import sys
import re
import socket
import copy
import logging

from pymatgen.io.vaspio.vasp_input import Incar, Poscar
from pymatgen.io.vaspio.vasp_input import Potcar, Kpoints
from pymatgen.apps.borg.queen import BorgQueen

from fireworks.core.firework import FireTaskBase, Firework, FWAction
from fireworks.core.launchpad import LaunchPad
from fireworks.utilities.fw_serializers import FWSerializable
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from mpinterfaces.calibrate import CalibrateMolecule, CalibrateSlab
from mpinterfaces.calibrate import CalibrateBulk
from mpinterfaces.measurement import Measurement
from mpinterfaces.database import MPINTVaspToDbTaskDrone

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)

def load_class(mod, name):
    """
    load class named name from module mod
    this function is adapted from materialsproject
    """
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)

def get_qadapter(nnodes=1, nprocs=16, walltime='24:00:00',
                 job_bin=None):
    d = {}
    #hipergator
    if 'ufhpc' in socket.gethostname():
        if job_bin is None:
            job_bin='/home/km468/Software/VASP/vasp.5.3.5/vasp'
        else:
            job_bin = job_bin
        d = {'type':'PBS',
             'params':
                 {
                'nnodes': str(nnodes),
                'ppnode': str(int(nprocs/nnodes)),
                'walltime': walltime,
                'job_name': 'vasp_job',
                'email': 'mpinterfaces@gmail.com',
                'notification_options': 'ae',
                'pre_rocket': '#PBS -l pmem=1000mb',
                'rocket_launch': 'mpirun '+job_bin
                }
             }
    #stampede
    elif 'stampede' in socket.gethostname():
        if job_bin is None:
            job_bin='/home1/01682/km468/Software/VASP/vasp.5.3.5/vasp'
        else:
            job_bin = job_bin
        d = {'type':'SLURM',
             'params':
                 {
                'nodes': str(nnodes),
                'ntasks': str(nprocs),
                'walltime': walltime,
                'queue':'normal',
                'account':'TG-DMR050028N',
                'job_name': 'vasp_job',
                'rocket_launch': 'ibrun '+job_bin
                }
             }
    if d:
        return CommonAdapter(d['type'], **d['params'])
    else:
        return None
    
def get_cal_obj(d):
    """
    construct a calibration object from the input dictionary, d

    returns a calibration object
    """
    incar = Incar.from_dict(d["incar"])
    poscar = Poscar.from_dict(d["poscar"])
    symbols = poscar.site_symbols
    potcar = Potcar(symbols)
    kpoints = Kpoints.from_dict(d["kpoints"])
    turn_knobs = d["turn_knobs"]
    #if running on hipergator or stampede, will return a qadapter
    #with default values even if que_params is not set
    qadapter = get_qadapter(**d.get('que_params', {}))
    job_cmd = None
    #running henniggroup machines
    if qadapter is None and  socket.gethostname() in ['hydrogen',
                                                      'helium',
                                                      'lithium',
                                                      'beryllium',
                                                      'carbon']:
        job_cmd = ['nohup', '/opt/openmpi_intel/bin/mpirun',
                   '-n', '24',
                   '/home/km468/Software/VASP/vasp.5.3.5/vasp']
        if d.get("other_params") is not None:
            d.get("other_params").update({'wait':False})
        else:
            d["other_params"] = {'wait':False}
    #if running on your desktop or laptop
    else:
        job_cmd=['ls', '-lt']
    if d.get('job_cmd'):
        job_cmd = d.get('job_cmd')
    cal =  load_class("mpinterfaces.calibrate",
                      d["calibrate"])(incar, poscar, potcar, kpoints,
                                      system = d.get("system", None),
                                      qadapter=qadapter, job_cmd=job_cmd,
                                      turn_knobs=turn_knobs,
                                      **d.get("other_params", {}))
    if d.get('job_dir_list'):
        cal.job_dir_list = d.get('job_dir_list')
    if d.get('job_ids'):
        cal.job_ids = d.get('job_ids')        
    return cal

@explicit_serialize
class MPINTCalibrateTask(FireTaskBase, FWSerializable):
    """
    Calibration Task
    """
    required_params = ["incar", "poscar", "kpoints", "calibrate",
                        "turn_knobs"]
    optional_params = ["que_params", "job_cmd", "system", "other_params"]

    def run_task(self, fw_spec):
        """
        launch jobs to the queue
        """
        cal = get_cal_obj(self)
        cal.setup()
        cal.run()
        d = cal.as_dict()
        d.update({'que_params':self.get('que_params')})
        return FWAction(mod_spec=[{'_push': {'cal_objs':d}}])
        
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
        logger.info('The measurement task will be constructed from {} calibration objects'
                    .format(len(fw_spec['cal_objs'])) )
        for calparams in fw_spec['cal_objs']:
            calparams.update({'que_params':self.get('que_params')})
            cal = get_cal_obj(calparams)
            cal_objs.append(cal)
        done = load_class("mpinterfaces.calibrate", "Calibrate").check_calcs(cal_objs)
        if not done:
            logger.info('Calibration not done yet. Try again later')
            logger.info('All subsequent fireworks will be defused')
            logger.info("""Try re-running this firework again later.
            Re-running this firework will activate all the subsequent foreworks too""")
            logger.info('This fireworks id = {}'.format(self.get("fw_id")))
            return FWAction(defuse_children=True)
            ### to enable dynamic workflow, uncomment the following
            #if self.get("fw_id"):
            #    fw_id = int(self.get("fw_id")) + 1
            #    self["fw_id"] = fw_id
            #    new_fw = Firework(MPINTMeasurementTask(self),
            #                      spec={'cal_objs':fw_spec['cal_objs']},
            #                      name = 'new_fw', fw_id = -fw_id)
            #else:
            #    new_fw = Firework(MPINTMeasurementTask(self),
            #                      spec={'cal_objs':fw_spec['cal_objs']},
            #                      name = 'new_fw')
            #    
            #return FWAction(detours=new_fw)
        else:
            measure = load_class("mpinterfaces.measurement",self['measurement'])(cal_objs,
                                                                                 **self.get("other_params", {}))
            job_cmd = None
            if self.get("job_cmd", None) is not None:
                job_cmd = self.get("job_cmd")
            measure.setup()
            measure.run(job_cmd = job_cmd)
            cal_list = []
            for cal in measure.cal_objs:
                d = cal.as_dict()
                d.update({'que_params':self.get('que_params')})
                cal_list.append(d)
            return FWAction(update_spec={'cal_objs':cal_list})

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
        drone = MPINTVaspToDbTaskDrone(**d.get("dbase_params", {}))
        queen = BorgQueen(drone)#, number_of_drones=ncpus)
        queen.serial_assimilate(self["measure_dir"])
        return FWAction()
