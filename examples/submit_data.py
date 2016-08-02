# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from pymatgen.apps.borg.queen import BorgQueen

from mpinterfaces.database import MPINTVaspToDbTaskDrone

# import multiprocessing

additional_fields = {"author": "kiran"}  # "doi":"10.1063/1.4865107"
drone = MPINTVaspToDbTaskDrone(host="127.0.0.1", port=27017,
                               database="vasp", collection="collection_name",
                               user="username", password="password",
                               additional_fields=additional_fields)

ncpus = 4  # multiprocessing.cpu_count()
queen = BorgQueen(drone, number_of_drones=ncpus)
queen.parallel_assimilate('path_to_vasp_calculation_folders')
