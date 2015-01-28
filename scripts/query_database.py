from __future__ import division, unicode_literals, print_function

from matgendb.query_engine import QueryEngine

qe = QueryEngine(host="127.0.0.1", port=27017,
                 database="vasp", collection="nanoparticles",
                 user="km468", password="km468")
results = list( qe.query(criteria={"hkl": [1,1,1]},
                         properties=['author','energy', 'ligand',
                                     'cif', 'pretty_formula']) )

print('formula : {}'.format(results[0]['pretty_formula']))
print('author : {}'.format(results[0]['author']))
print('ligand : {}'.format(results[0]['ligand']))
print('energy : {}'.format(results[0]['energy']))
print('CIF structure data = \n {}'.format(results[0]['cif']))

