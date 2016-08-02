# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os

from matgendb.query_engine import QueryEngine
from pymatgen.core.structure import Structure

# from config file db.json
DB_CONFIG = os.path.join(os.path.expanduser('~'), ".mongodb/db.json")
qe = QueryEngine.from_config(DB_CONFIG)
# or 
# qe = QueryEngine(host="127.0.0.1", port=27017,
#                 database="vasp", collection="collection_name",
#                 user="username", password="password")

results = qe.query(criteria={"normalized_formula": 'GaSb'},
                   properties=['pretty_formula', 'author', 'spacegroup',
                               'output', 'analysis', 'last_updated',
                               'dir_name'])
# if the documents contain the hkl field, to query based on hkl
# use criteria={"hkl": [1,1,1]}

for r in results:
    for k, v in r.items():
        if k == "output":
            structure = Structure.from_dict(v["crystal"])
            print(structure)
            print('\nFinal energy : {}\n'.format(v["final_energy"]))
        elif k == "analysis":
            print('Band gap : {}\n'.format(v["bandgap"]))
        else:
            print('{0} : \n{1}\n'.format(k, v))
