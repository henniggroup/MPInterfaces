# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
process vasprun.xml file by walking through the enitre directory tree
in the parent directory
"""

import os
import glob
import re

from monty.json import MontyDecoder

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.entries.computed_entries import ComputedEntry, \
    ComputedStructureEntry
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.hive import _get_transformation_history

from mpinterfaces.default_logger import get_default_logger

logger = get_default_logger(__name__)


class MPINTComputedEntry(ComputedEntry):
    """

    extend ComputedEntry to include structure as well as kpoints

    """

    def __init__(self, structure, kpoints, incar, energy, correction=0.0,
                 parameters=None, data=None, entry_id=None):
        ComputedEntry.__init__(self, structure.composition, energy,
                               correction=correction,
                               parameters=parameters,
                               data=data, entry_id=entry_id)
        self.structure = structure
        self.kpoints = kpoints
        self.incar = incar
        # self.data = {"style": self.kpoints.style,
        # 'kpoints': self.kpoints.kpts, 'incar':self.incar.as_dict()}

    def __repr__(self):
        output = ["MPINTComputedEntry {}".format(self.composition.formula),
                  "Energy = {:.4f}".format(self.uncorrected_energy),
                  "Correction = {:.4f}".format(self.correction), "Parameters:"]
        for k, v in self.parameters.items():
            output.append("{} = {}".format(k, v))
        output.append("Data:")
        for k, v in self.data.items():
            output.append("{} = {}".format(k, v))
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()

    def as_dict(self):
        d = super(ComputedEntry, self).as_dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["structure"] = self.structure.as_dict()
        d["kpoints"] = self.kpoints.as_dict()
        d["incar"] = self.incar.as_dict()
        return d

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        return cls(dec.process_decoded(d["structure"],
                                       d["kpoints"],
                                       d["incar"]),
                   d["energy"], d["correction"],
                   dec.process_decoded(d.get("parameters", {})),
                   dec.process_decoded(d.get("data", {})),
                   entry_id=d.get("entry_id", None))


class MPINTVasprun(Vasprun):
    """
    Extend Vasprun to use custom ComputedEntry: MPINTComputedEntry
    """

    def __init__(self, filename, ionic_step_skip=None,
                 ionic_step_offset=0, parse_dos=True,
                 parse_eigen=True, parse_projected_eigen=False,
                 parse_potcar_file=True):

        Vasprun.__init__(self, filename,
                         ionic_step_skip=ionic_step_skip,
                         ionic_step_offset=ionic_step_offset,
                         parse_dos=parse_dos, parse_eigen=parse_eigen,
                         parse_projected_eigen=parse_projected_eigen,
                         parse_potcar_file=parse_potcar_file)

    def get_computed_entry(self, inc_structure=False,
                           inc_incar_n_kpoints=False,
                           parameters=None, data=None):
        """
        Returns a ComputedEntry from the vasprun.

        Args:
            inc_structure (bool): Set to True if you want
                ComputedStructureEntries to be returned instead of
                ComputedEntries.
            inc_incar_n_kpoints (bool): along with inc_structure set
                to True if you want MPINTComputedEntries to be
                returned
            parameters (list): Input parameters to include.
                It has to be one of the properties supported by the
                Vasprun object. If parameters == None, a default set
                of parameters that are necessary for typical
                post-processing will be set.
            data (list): Output data to include. Has to be one of the
                properties supported by the Vasprun object.

        Returns:
            ComputedStructureEntry/ComputedEntry
        """
        param_names = {"is_hubbard", "hubbards", "potcar_symbols",
                       "run_type"}
        if parameters:
            param_names.update(parameters)
        params = {p: getattr(self, p) for p in param_names}
        data = {p: getattr(self, p) for p in data} if data is not None else {}

        if inc_structure and inc_incar_n_kpoints:
            return MPINTComputedEntry(self.final_structure,
                                      self.kpoints, self.incar,
                                      self.final_energy,
                                      parameters=params, data=data)

        elif inc_structure:
            return ComputedStructureEntry(self.final_structure,
                                          self.final_energy,
                                          parameters=params,
                                          data=data)
        else:
            return ComputedEntry(self.final_structure.composition,
                                 self.final_energy, parameters=params,
                                 data=data)


class MPINTVaspDrone(VaspToComputedEntryDrone):
    """
    extend VaspToComputedEntryDrone to use the custom
    Vasprun: MPINTVasprun
    """

    def __init__(self, inc_structure=False, inc_incar_n_kpoints=False,
                 parameters=None, data=None):
        VaspToComputedEntryDrone.__init__(self,
                                          inc_structure=inc_structure,
                                          parameters=parameters,
                                          data=data)
        self._inc_structure = inc_structure
        self._inc_incar_n_kpoints = inc_incar_n_kpoints
        self._parameters = parameters
        self._data = data

    def assimilate(self, path):
        files = os.listdir(path)
        if "relax1" in files and "relax2" in files:
            filepath = glob.glob(os.path.join(path, "relax2",
                                              "vasprun.xml*"))[0]
        else:
            vasprun_files = glob.glob(os.path.join(path, "vasprun.xml*"))
            filepath = None
            if len(vasprun_files) == 1:
                filepath = vasprun_files[0]
            elif len(vasprun_files) > 1:
                for fname in vasprun_files:
                    if os.path.basename(fname) in ["vasprun.xml",
                                                   "vasprun.xml.gz",
                                                   "vasprun.xml.bz2"]:
                        filepath = fname
                        break
                    if re.search("relax2", fname):
                        filepath = fname
                        break
                    filepath = fname
        try:
            vasprun = MPINTVasprun(filepath)
        except Exception as ex:
            logger.error("vasprun read error in {}: {}".format(filepath, ex))
            return None
        entry = vasprun.get_computed_entry(self._inc_structure,
                                           self._inc_incar_n_kpoints,
                                           parameters=self._parameters,
                                           data=self._data)
        entry.parameters["history"] = _get_transformation_history(path)
        return entry

    def __str__(self):
        return " MPINTVaspDrone"

    def as_dict(self):
        return {"init_args": {"inc_structure": self._inc_structure,
                              "inc_incar_n_kpoints": self._inc_incar_n_kpoints,
                              "parameters": self._parameters,
                              "data": self._data},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])
