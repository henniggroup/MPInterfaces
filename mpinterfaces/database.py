from __future__ import division, unicode_literals, print_function

"""
Put data into mongo database
"""

import sys
import os
import json
import logging

from pymatgen.io.vaspio import Outcar

from matgendb.creator import VaspToDbTaskDrone, get_uri

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


class MPINTVaspToDbTaskDrone(VaspToDbTaskDrone):  
    """
    subclassing VaspToDbTaskDrone
    """
    def __init__(self, host="10.5.46.101", port=27017, database="vasp",
                 user=None, password=None, collection="nanoparticles",
                 parse_dos=False, compress_dos=False, simulate_mode=False,
                 additional_fields=None, update_duplicates=True,
                 mapi_key=None, use_full_uri=True, runs=None):
        """
        Additional_fields = {'authors':[a.as_dict() for a in self.authors],
                             'projects': ,
                              'refernces': ,
                              }
        """
        VaspToDbTaskDrone.__init__(self, host=host, port=port, 
                                   database=database, user=user, 
                                   password=password, 
                                   collection=collection,
                                   parse_dos=parse_dos, 
                                   compress_dos=compress_dos,
                                   simulate_mode=simulate_mode, 
                                   additional_fields=additional_fields,
                                   update_duplicates=update_duplicates,
                                   mapi_key=mapi_key, 
                                   use_full_uri=use_full_uri, runs=runs)


    def post_process(self, dir_name, d):
        """
        customization:
            adds system.json to the dictionary
        """
        logger.info("Post-processing dir:{}".format(dir_name))
        fullpath = os.path.abspath(dir_name)
        filename = os.path.join(fullpath, "system.json")
        with open(filename, "r") as f:
            system = json.load(f)
            d["hkl"] = system["hkl"]
            d["ligand"] = system["ligand"]
        #try:
        #    run_stats = {}
        #    outcar = Outcar("OUTCAR")
        #    #print(outcar.efermi)
        #    taskname = "run1"
        #    d["calculations"][0]["output"]["outcar"] = outcar.as_dict()
        #    run_stats[taskname] = outcar.run_stats
        #except:
        #    logger.error("Bad OUTCAR for {}.".format(fullpath))
        #try:
        #    overall_run_stats = {}
        #    for key in ["Total CPU time used (sec)", "User time (sec)",
        #                "System time (sec)", "Elapsed time (sec)"]:
        #        overall_run_stats[key] = sum([v[key]
        #                                      for v in run_stats.values()])
        #    run_stats["overall"] = overall_run_stats
        #except:
        #    logger.error("Bad run stats for {}.".format(fullpath))
        #d["run_stats"] = run_stats
        #Convert to full uri path.
        if self.use_full_uri:
            d["dir_name"] = get_uri(dir_name)
        logger.info("Post-processed " + fullpath)
