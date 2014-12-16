#4 fireworks: 1 firetask each

#- calibration: fwaction return list of dirs
#- get knob settings: fwaction returns the dict of incar, poscar, potcar, kpoints objects

#- run experiments: fwaction returns list of run dirs
#- to database: make the measurements(binding enrgy ect) and put it in the database
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.core.launchpad import LaunchPad
from fireworks.utilities.fw_serializers import FWSerializable

from matgendb.creator import VaspToDbTaskDrone


class MPINTCalibrateTask(FireTaskBase, FWSerializable):
    pass



class MPINTMeasurementRunTask(FireTaskBase, FWSerializable):
    pass



class MPINTMeasurementTask(FireTaskBase, FWSerializable):
    pass



class MPINTVaspToDBTask(FireTaskBase, FWSerializable):

    _fw_name = "MPINT Vasp to Database Task"


    def run_task(self, fw_spec):
        prev_dir = fw_spec['prev_vasp_dir']

        # get the db credentials                                                                                            
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'tasks_db.json')

        # use MPDrone to put it in the database                                                                             
        with open(db_path) as f:
            db_creds = json.load(f)
            drone = VaspToDbTaskDrone(
                host=db_creds['host'], port=db_creds['port'],
                database=db_creds['database'], user=db_creds['admin_user'],
                password=db_creds['admin_password'],
                collection=db_creds['collection'])
            t_id = drone.assimilate('Measurement')

        if t_id:
            print 'ENTERED task id:', t_id
            stored_data = {'task_id': t_id}
            update_spec = {'prev_vasp_dir': prev_dir, 'prev_task_type': fw_spec['prev_task_type']}
            return FWAction(stored_data=stored_data, update_spec=update_spec)
        else:
            raise ValueError("Could not parse entry for database insertion!")
