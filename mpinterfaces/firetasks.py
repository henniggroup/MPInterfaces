"""
Defines various firetasks
"""

import re
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from mpinterfaces.calibrate import CalibrateMolecule
from mpinterfaces.measurement import Measurement
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.core.launchpad import LaunchPad
from fireworks.utilities.fw_serializers import FWSerializable
from fireworks.utilities.fw_utilities import explicit_serialize
from matgendb.creator import VaspToDbTaskDrone
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter


def load_class(mod, name):
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)


def to_int_list(str_list):
    m = re.search(r"\[(\d+)\,(\d+)\,(\d+)\]", str_list)
    return [int(m.group(i)) for i in range(1,4)]       


def get_cal_obj(d):
    qadapter = None
    if d['que']:
        qadapter = CommonAdapter(d['que']['type'], **d['que']['params'])
    incar = Incar.from_dict(d["incar"])
    poscar = Poscar.from_dict(d["poscar"])
    symbols = poscar.site_symbols #symbol_set
    potcar = Potcar(symbols)
    kpoints = Kpoints.from_dict(d["kpoints"])
    if qadapter is not None:
        cal =  load_class("mpinterfaces.calibrate", d["calibrate"])(incar, poscar, potcar, kpoints, qadapter=qadapter, job_cmd='qsub', **d.get("cal_construct_params", {}))
    elif d['job_cmd']:
        cal =  load_class("mpinterfaces.calibrate", d["calibrate"])(incar, poscar, potcar, kpoints, job_cmd=d['job_cmd'], **d.get("cal_construct_params", {}))
    #no qadapter and no job_cmd
    else:
        cal =  load_class("mpinterfaces.calibrate", d["calibrate"])(incar, poscar, potcar, kpoints, job_cmd='vasp', **d.get("cal_construct_params", {}))            
    return cal           
    


@explicit_serialize
class MPINTCalibrateTask(FireTaskBase, FWSerializable):
    """
    incar: Incar object.to_dict
    similarly for poscar
    encut_list : example:- ['400', '800', '100'] --> range specification
    kpoint_list: example:- ['[7,7,7]', '[11,11,11]' ]
    """
    
    required_params = ["incar", "poscar", "kpoints","calibrate", "que"]
    optional_params = ["job_cmd", "encut_list", "kpoint_list", "cal_construct_params"]
#    _fw_name = 'MPINTCalibrateTask'

    def run_task(self, fw_spec):
        """
        launch jobs to the queue
        """
        cal = get_cal_obj(self)
        if self["encut_list"]:
            range_specs = [ int(encut) for encut in self["encut_list"] ]
            encut_list = range(range_specs[0], range_specs[1], range_specs[2])
            cal.encut_cnvg(encut_list)
        if self["kpoint_list"]:
            kpoint_list = [ to_int_list(kpt) for kpt in self["kpoint_list"]]
            cal.kpoints_cnvg(kpoints_list=kpoint_list)        
        cal.run()
        

@explicit_serialize
class MPINTMeasurementTask(FireTaskBase, FWSerializable):
    
    required_params = ["cal_objs"]
    optional_params = ["msr_construct_params"]    
#    _fw_name = 'MPINTMeasurementTask'

    def run_task(self, fw_spec):
        """
        go through the calibration directiories and get the optimu knob_settings
        and launch the actual measurement jobs to the queue
        """
        cal_objs_list = []
        for calparams in self['cal_objs']:
            cal = get_cal_obj(calparams)
            cal_objs_list.append(cal)
        measure = Measurement(cal_objs_list, **self.get("msr_construct_params", {}))
        #test
        #measure.setup()
        measure.calmol.knob_settings('1')


@explicit_serialize
class MPINTPostProcessTask(FireTaskBase, FWSerializable):

    required_params = ["measure_dirs"]    
#    _fw_name = "MPINTPostProcessTask"


    def run_task(self, fw_spec):
        """
        go through the measurement job dirs and compute quatities of interst
        also put the measurement jobs in the datbase
        """

        drone = VaspToDbTaskDrone(
                host='localhost', port='27017',
                database='vasp_test', user='km468',
                password='km468',
                collection='test_collection')
        t_id = drone.assimilate('Measurement')

        if t_id:
            print 'ENTERED task id:', t_id
            stored_data = {'task_id': t_id}
            update_spec = {'prev_vasp_dir': prev_dir, 'prev_task_type': fw_spec['prev_task_type']}
            return FWAction(stored_data=stored_data, update_spec=update_spec)
        else:
            raise ValueError("Could not parse entry for database insertion!")


#        # get the db credentials                                                             
#        db_dir = os.environ['DB_LOC']
#        db_path = os.path.join(db_dir, 'tasks_db.json')
#
#        # use MPDrone to put it in the database                                              #                               
#        with open(db_path) as f:
#            db_creds = json.load(f)
#            drone = VaspToDbTaskDrone(
#                host=db_creds['host'], port=db_creds['port'],
#                database=db_creds['database'], user=db_creds['admin_user'],
#                password=db_creds['admin_password'],
#                collection=db_creds['collection'])
#            t_id = drone.assimilate('Measurement')
        
