from __future__ import division, unicode_literals, print_function

from mpinterfaces.database import MPINTVaspToDbTaskDrone
from pymatgen.apps.borg.queen import BorgQueen
#import multiprocessing

additional_fields = {"author":"kiran"}
drone = MPINTVaspToDbTaskDrone(host="127.0.0.1", port=27017,
                               database="vasp", collection="collection_name",
                               user="username", password="password",
                               additional_fields=additional_fields)

ncpus = 4 #multiprocessing.cpu_count()
queen = BorgQueen(drone, number_of_drones=ncpus)
queen.parallel_assimilate('path_to_vasp_calculation_folders')
