from __future__ import division, unicode_literals, print_function

import os
from matgendb.query_engine import QueryEngine
from monty.json import MontyDecoder
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# from config file db.json
DB_CONFIG = os.path.join(os.path.expanduser('~'), ".mongodb/db.json")
qe = QueryEngine.from_config(DB_CONFIG)
# or 
#qe = QueryEngine(host="127.0.0.1", port=27017,
#                 database="vasp", collection="collection_name",
#                 user="username", password="password")

results = qe.query(criteria = {"pretty_formula": 'InSb'},
                   properties = ['pretty_formula','author',
                                 'input', 'output'])
# if the documents contain the hkl field, to query based on hkl
# use criteria={"hkl": [1,1,1]}

for r in results:
    for k, v in r.items():
        print('{0} : \n{1}\n'.format(k,v))
# convert to pymatgen structure and get the spacegroup
        if k == "output":
            structure = MontyDecoder().process_decoded(v["crystal"])
            sga = SpacegroupAnalyzer(structure)
            print("Final structure space group: {}".format(sga.get_spacegroup_symbol()))
