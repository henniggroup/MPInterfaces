# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Put data into mongo database
"""

from six.moves import range

import sys
import os
import json
import logging
import socket
import string
import datetime

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.bond_valence import BVAnalyzer

from matgendb.creator import VaspToDbTaskDrone
from matgendb.creator import logger as mgdb_logger

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)
mgdb_logger.addHandler(sh)


class MPINTVaspToDbTaskDrone(VaspToDbTaskDrone):
    """
    subclassing VaspToDbTaskDrone
    """

    def __init__(self, host="127.0.0.1", port=27017, database="vasp",
                 user=None, password=None, collection="nanoparticles",
                 parse_dos=False, compress_dos=False,
                 simulate_mode=False,
                 additional_fields=None, update_duplicates=True,
                 mapi_key=None, use_full_uri=True, runs=None):
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
                                   use_full_uri=use_full_uri,
                                   runs=runs)

    def generate_doc(self, dir_name, vasprun_files):
        """
        Overridden
        """
        try:
            fullpath = os.path.abspath(dir_name)
            d = {k: v for k, v in self.additional_fields.items()}
            d["dir_name"] = fullpath
            d["schema_version"] = VaspToDbTaskDrone.__version__
            d["calculations"] = [
                self.process_vasprun(dir_name, taskname, filename)
                for taskname, filename in vasprun_files.items()]
            d1 = d["calculations"][0]
            d2 = d["calculations"][-1]
            # Now map some useful info to the root level.
            for root_key in ["completed_at", "nsites",
                             "unit_cell_formula",
                             "reduced_cell_formula",
                             "pretty_formula",
                             "elements", "nelements", "cif",
                             "density",
                             "is_hubbard", "hubbards", "run_type"]:
                d[root_key] = d2[root_key]
            d["chemsys"] = "-".join(sorted(d2["elements"]))
            # store any overrides to the exchange correlation functional
            xc = d2["input"]["incar"].get("GGA")
            if xc:
                xc = xc.upper()
            d["input"] = {"crystal": d1["input"]["crystal"],
                          "is_lasph": d2["input"]["incar"].get("LASPH", False),
                          "potcar_spec": d1["input"].get("potcar_spec"),
                          "xc_override": xc}
            vals = sorted(d2["reduced_cell_formula"].values())
            d["anonymous_formula"] = {string.ascii_uppercase[i]: float(vals[i])
                                      for i in range(len(vals))}
            d["output"] = {
                "crystal": d2["output"]["crystal"],
                "final_energy": d2["output"]["final_energy"],
                "final_energy_per_atom": d2["output"]["final_energy_per_atom"]}
            d["name"] = "vasp"
            p = d2["input"]["potcar_type"][0].split("_")
            pot_type = p[0]
            functional = "lda" if len(pot_type) == 1 else "_".join(p[1:])
            d["pseudo_potential"] = {"functional": functional.lower(),
                                     "pot_type": pot_type.lower(),
                                     "labels": d2["input"]["potcar"]}
            if len(d["calculations"]) == len(self.runs) or \
                    list(vasprun_files.keys())[0] != "relax1":
                d["state"] = "successful" if d2["has_vasp_completed"] \
                    else "unsuccessful"
            else:
                d["state"] = "stopped"
            d["analysis"] = analysis_and_error_checks(d)
            sg = SpacegroupAnalyzer(
                Structure.from_dict(d["output"]["crystal"]), 0.1)
            d["spacegroup"] = {"symbol": sg.get_space_group_symbol(),
                               "number": sg.get_space_group_number(),
                               "point_group": sg.get_point_group(),
                               "source": "spglib",
                               "crystal_system": sg.get_crystal_system(),
                               "hall": sg.get_hall()}
            d["last_updated"] = datetime.datetime.today()
            return d
        except Exception as ex:
            import traceback
            print(traceback.format_exc())
            logger.error("Error in " + os.path.abspath(dir_name) +
                         ".\n" + traceback.format_exc())
            return None

    def post_process(self, dir_name, d):
        """
        customization:
            adds system.json to the dictionary
        """
        logger.info("Post-processing dir:{}".format(dir_name))
        fullpath = os.path.abspath(dir_name)
        filename = os.path.join(fullpath, "system.json")
        if os.path.exists(filename):
            with open(filename, "r") as f:
                system = json.load(f)
                d["hkl"] = system.get("hkl")
                d["ligand"] = system.get("ligand")
        # from pyamtgen-db
        # Parse OUTCAR for additional information and run
        # stats that are generally not in vasprun.xml.
        try:
            run_stats = {}
            overall_run_stats = {}
            for key in ["Total CPU time used (sec)", "User time (sec)",
                        "System time (sec)", "Elapsed time (sec)"]:
                overall_run_stats[key] = sum([v[key]
                                              for v in run_stats.values()])
            run_stats["overall"] = overall_run_stats
        except:
            logger.error("Bad run stats for {}.".format(fullpath))
        d["run_stats"] = run_stats
        # Convert to full uri path.
        if self.use_full_uri:
            d["dir_name"] = get_uri(dir_name)
        logger.info("Post-processed " + fullpath)


def get_uri(dir_name):
    """
    Customized version of the original pymatgen-db version.
    Customization required because same job folder on hipergator
    gets different uri for different login nodes .

    Returns the URI path for a directory. This allows files hosted on
    different file servers to have distinct locations.
    Args:
        dir_name:
            A directory name.
    Returns:
        Full URI path, e.g., fileserver.host.com:/full/path/of/dir_name.
    """
    fullpath = os.path.abspath(dir_name)
    try:
        hostname = socket.gethostbyaddr(socket.gethostname())[0].split('.')[1]
    except:
        hostname = socket.gethostname()
    return "{}:{}".format(hostname, fullpath)


# remove coordination number
def analysis_and_error_checks(d, max_force_threshold=0.5,
                              volume_change_threshold=0.2):
    initial_vol = d["input"]["crystal"]["lattice"]["volume"]
    final_vol = d["output"]["crystal"]["lattice"]["volume"]
    delta_vol = final_vol - initial_vol
    percent_delta_vol = delta_vol / initial_vol
    # coord_num = get_coordination_numbers(d)
    calc = d["calculations"][-1]
    gap = calc["output"]["bandgap"]
    cbm = calc["output"]["cbm"]
    vbm = calc["output"]["vbm"]
    is_direct = calc["output"]["is_gap_direct"]

    warning_msgs = []
    error_msgs = []

    if abs(percent_delta_vol) > volume_change_threshold:
        warning_msgs.append("Volume change > {}%"
                            .format(volume_change_threshold * 100))

    bv_struct = Structure.from_dict(d["output"]["crystal"])
    try:
        bva = BVAnalyzer()
        bv_struct = bva.get_oxi_state_decorated_structure(bv_struct)
    except ValueError as e:
        logger.error("Valence cannot be determined due to {e}."
                     .format(e=e))
    except Exception as ex:
        logger.error("BVAnalyzer error {e}.".format(e=str(ex)))

    max_force = None
    if d["state"] == "successful" and \
        d["calculations"][0]["input"]["parameters"].get("NSW",
                                                        0) > 0:
        # handle the max force and max force error
        max_force = max([np.linalg.norm(a)
                         for a in d["calculations"][-1]["output"]
                         ["ionic_steps"][-1]["forces"]])

        if max_force > max_force_threshold:
            error_msgs.append("Final max force exceeds {} eV"
                              .format(max_force_threshold))
            d["state"] = "error"

        s = Structure.from_dict(d["output"]["crystal"])
        if not s.is_valid():
            error_msgs.append("Bad structure (atoms are too close!)")
            d["state"] = "error"

    return {"delta_volume": delta_vol,
            "max_force": max_force,
            "percent_delta_volume": percent_delta_vol,
            "warnings": warning_msgs,
            "errors": error_msgs,
            "bandgap": gap, "cbm": cbm, "vbm": vbm,
            "is_gap_direct": is_direct,
            "bv_structure": bv_struct.as_dict()}
