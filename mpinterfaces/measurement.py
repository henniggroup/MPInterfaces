"""

combines instrument, calibrate and interfaces to perform the calibration and run the actual jobs

"""

import os
import operator
import glob
from collections import Counter
import numpy as np
from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vaspio.vasp_output import Vasprun
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone, SimpleVaspToComputedEntryDrone, _get_transformation_history
from pymatgen.apps.borg.queen import BorgQueen
#pymatgen.entries.computed_entries
from mpinterfaces.calibrate import Calibrate
from mpinterfaces.interface import Interface


class MPINTComputedEntry(ComputedEntry):
        
    """
        
    extend ComputedEntry to include structure as well as kpoints
    
    """

    def __init__(self, structure, kpoints, incar, energy, correction=0.0,
                  parameters=None, data=None, entry_id=None):
        ComputedEntry.__init__(self, structure.composition, energy,
                               correction=correction, parameters=parameters,
                               data=data, entry_id=entry_id)
        self.structure = structure
        self.kpoints = kpoints
        self.incar = incar
        #self.data = {"style": self.kpoints.style, 'kpoints': self.kpoints.kpts, 'incar':self.incar.as_dict()}
        

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
        return cls(dec.process_decoded(d["structure"],d["kpoints"],d["incar"]),
                   d["energy"], d["correction"],
                   dec.process_decoded(d.get("parameters", {})),
                   dec.process_decoded(d.get("data", self.data)),
                   entry_id=d.get("entry_id", None))




class MPINTVasprun(Vasprun):
        
    """
        
    Extend Vasprun to use custom ComputedEntry: MPINTComputedEntry
    
    """
        
    def __init__(self, filename, ionic_step_skip=None,
                 ionic_step_offset=0, parse_dos=True,
                 parse_eigen=True, parse_projected_eigen=False):
        
        Vasprun.__init__(self, filename, ionic_step_skip=ionic_step_skip,
                 ionic_step_offset=ionic_step_offset, parse_dos=parse_dos,
                 parse_eigen=parse_eigen, parse_projected_eigen=parse_projected_eigen)
            

    def get_computed_entry(self, inc_structure=False, inc_incar_n_kpoints=False,
                           parameters=None, data=None):
        """
        Returns a ComputedEntry from the vasprun.

        Args:
            inc_structure (bool): Set to True if you want
                ComputedStructureEntries to be returned instead of
                ComputedEntries.
            inc_incar_n_kpoints (bool): along with inc_structure set to True if you want
                MPINTComputedEntries to be returned : returns incar and kpoints objects 
                
            parameters (list): Input parameters to include. It has to be one of
                the properties supported by the Vasprun object. If
                parameters == None, a default set of parameters that are
                necessary for typical post-processing will be set.
            data (list): Output data to include. Has to be one of the properties
                supported by the Vasprun object.

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
            return MPINTComputedEntry(self.final_structure, self.kpoints, self.incar, 
                                        self.final_energy, parameters=params,
                                          data=data)
        
        elif inc_structure:
            return ComputedStructureEntry(self.final_structure,
                                          self.final_energy, parameters=params,
                                          data=data)
        else:
            return ComputedEntry(self.final_structure.composition,
                                 self.final_energy, parameters=params,
                                 data=data)



class MPINTVaspDrone(VaspToComputedEntryDrone):
        
    """
    
    extend VaspToComputedEntryDrone to use the custom Vasprun: MPINTVasprun
        
    """

    def __init__(self, inc_structure=False,  inc_incar_n_kpoints=False, parameters=None, data=None):
        VaspToComputedEntryDrone.__init__(self, inc_structure=inc_structure, parameters=parameters, data=data)
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
            print "error in {}: {}".format(filepath, ex)
            #logger.debug("error in {}: {}".format(filepath, ex))
            return None

        entry = vasprun.get_computed_entry(self._inc_structure, self._inc_incar_n_kpoints,
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
                "version": __version__,
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])
    
    

        
        
        

class Measurement(Calibrate):
    
    """
    
    Perfor the actual measurement using the optimum knob settings
    takes calibratemolecule, calibrateslab and interface objects and
    perform the interface calculations and compute binding energy, bandstructure etc. 
    
    """

    def __init__(self, incar, poscar, potcar, kpoints, setup_dir='.', parent_job_dir='./Measurement'):
        Calibrate.__init__(self, incar, poscar, potcar, kpoints, setup_dir=setup_dir, parent_job_dir=parent_job_dir)
        #self.interface = interface
        self.encut = None
        self.kpoints = None
        self.vac_spacing = None
        self.slab_thickness = None
        self.jobs = []
        self.handlers = []


    def knob_settings(self, rootpath):
        """
        go through the parent dir and get all encut, kpoints and energies
        also vac spacing and slab thinckness for slab calulations
        these values willl be used to do the actual interface measurements
        use Vasprun class to get the afore mentioned values from the xml files
        should not proceed if the calculations are not done
        should update the incar, poscar, potcar, kpoints objects according to
        the knob_settings
        
        
        """
        encut = []
        kpts = []
        drone = MPINTVaspDrone(inc_structure=True, inc_incar_n_kpoints=True) #VaspToComputedEntryDrone()#
        bg =  BorgQueen(drone)
        bg.serial_assimilate(rootpath)
        alldata =  bg.get_data()
        enkp=[]
        c = Counter()
        for d in alldata:
            if d:
                enkp.append(str(d.incar['ENCUT']))
                enkp.append(str(d.kpoints.kpts))
                c[str(d.incar['ENCUT'])] += 1
                #encut.append( [ str(d.incar['ENCUT']), d.energy] )
                #kpts.append( [ str(d.kpoints.kpts),  d.energy] )
        enkp_mc =  Counter(enkp).most_common(2)
        kp_mc = None
        en_mc = None
        if '[[' in enkp_mc[0][0]:
            kp_mc = enkp_mc[0][0]
            en_mc = enkp_mc[1][0]
        else:
            kp_mc = enkp_mc[1][0]
            en_mc = enkp_mc[0][0]
        energy_kpt={}
        energy_encut = {}
        for d in alldata:
            if d:
                if str(d.incar['ENCUT']) == en_mc:
                    energy_kpt[str(d.kpoints.kpts)] = d.energy
                if str(str(d.kpoints.kpts)) == kp_mc:
                    energy_encut[str(d.incar['ENCUT'])] = d.energy
        print energy_encut
        print energy_kpt
        #from large to small
        sorted_encut = sorted(energy_encut.items(), key=operator.itemgetter(1), reverse=True)
        sorted_kpt = sorted(energy_kpt.items(), key=operator.itemgetter(1), reverse=True)
        print 'sorted'
        print sorted_encut
        print sorted_kpt
        delta_e = 0.01
        possible_encuts = []
        for i, e in enumerate(sorted_encut):
            print i, e
            if i < len(sorted_encut)-1:
                if np.abs(sorted_encut[i+1][1] - e[1]) < 0.01:
                    possible_encuts.append(sorted_encut[i+1][0])
        if possible_encuts:
            print possible_encuts
        else:
            print 'convergence not reached'


                    
        
       



    #def get_


    def setup(self, interfaces):
        """
        create a list of jobs to run
        exampls:- change coverage, site occupancies,
        interfaces: list of interfaces for which jobs must be created
        also run molecule, slab relaxation calcs
        """
        for iface in interfaces:
            job_dir  = self.parent_job_dir + iface.name
            self.job_dirs.append(job_dir)
            vis = MPINTVaspInputSet(iface.name, self.incar, self.poscar, self.potcar, self.kpoints)
            #the job command can be overrridden in the run method
            job = MPINTVaspJob(["pwd"], final = True, setup_dir=self.setup_dir, job_dir=job_dir, vis=vis, auto_npar=False, auto_gamma=False)
            self.jobs.append(job)

    def make_measurements(self):
        """
        To calculate binding energy of interface.
        Needs to get relaxed molecule energy, relaxed slab energy and relaxed interface energies
        Get the molecule and slab energies from a relaxation run
        of optimum paramter set for Molecule and Slab (when is
        that relaxation run done? NOTE: static done with
        as of Calibration, do before this call)
        Use Vasprun() to get energies, same method as knob_settings in Measurement.

        Goes to the required directories for Molecule and Slab and uses Vasprun to parse through the vasprun.xml and gets the energy, similar for 
        """
        
        pass
        



#test
if __name__=='__main__':
    system = 'Pt bulk'
    atoms = ['Pt']
    a0 = 3.965
    lvec = [ [0.5, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5] ]
    lvec = np.array(lvec) * a0
    lattice = Lattice(lvec)
    structure = Structure( lattice, atoms, [ [0.0, 0.0, 0.0] ],
                           coords_are_cartesian=False,
                           site_properties={"magmom":[0]} )

    incarparams = {'System':'test',
                   'ENCUT': 400,
                   'ISMEAR': 1,
                   'SIGMA': 0.1,
                   'EDIFF':1E-6}
    incar = Incar(params=incarparams)
    poscar = Poscar(structure, comment=system,
                    selective_dynamics=None,
                    true_names=True, velocities=None, predictor_corrector=None)
    potcar = Potcar(symbols=poscar.site_symbols, functional='PBE', sym_potcar_map=None)
    kpoints = Kpoints.monkhorst_automatic(kpts=(16, 16, 16), shift=(0, 0, 0))

    measure = Measurement(incar, poscar, potcar, kpoints)
    #get all data in all the directories in the provided rootfolder, here 1/
    measure.knob_settings('1')


    
    #print type(alldata)
    #for d in  alldata:
    #    if d:
    #        print 'final energy with out entropy ', d.energy
    #        #print 'kpoints : ', d.data['kpoints']
    #        #print 'incar : ', d.incar
    #        print 'kpoints = ', d.kpoints.kpts            
    #        print 'ENCUT = ', d.incar['ENCUT']                        
       
    
    
#    interfaces = Interface objects    
#    measure = Measurement(incar, poscar, potcar, kpoints)
#    measure.knob_settings()
#    measure.setup(interfaces)
#    measure.run(['qsub','job_script'])
