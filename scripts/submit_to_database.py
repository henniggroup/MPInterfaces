from __future__ import division, unicode_literals, print_function

from mpinterfaces.database import MPINTVaspToDbTaskDrone
from pymatgen.apps.borg.queen import BorgQueen
#import multiprocessing

additional_fields = {"author":"kiran"}
drone = MPINTVaspToDbTaskDrone(host="127.0.0.1", port=27017,
                               database="vasp", collection="nanoparticles",
                               user="km468", password="km468",
                               additional_fields=additional_fields)
#ncpus = multiprocessing.cpu_count()
#print ncpus
queen = BorgQueen(drone)#, number_of_drones=ncpus)
queen.serial_assimilate('/home/matk/Software/vasp_automation/test/Li2O')
