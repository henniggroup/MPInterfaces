from __future__ import division, unicode_literals, print_function

from matgendb.query_engine import QueryEngine

qe = QueryEngine(host="127.0.0.1", port=27017,
                 database="vasp", collection="collection_name",
                 user="username", password="password")

results = list( qe.query(criteria={"pretty_formula": 'InSb'},
                         properties=['input', 'output']) )

print('input : {}'.format(results[0]['input']))
print('output : {}'.format(results[0]['output']))

# if the documents contain the hkl field, to query based on hkl
# use criteria={"hkl": [1,1,1]}
