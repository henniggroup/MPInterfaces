from __future__ import division, unicode_literals, print_function

"""

ligand - slab interface calibration, relaxation, solvation workflows

1. get the initial structures from the materials project database
2. set relaxation calibration task for the hkl slab --> firework1
3. set static measurement task for each of the structures --> firework2
4. construct workflow consiting of the above mentioned fireworks for
   each of the structures and hkls
5. submit all the workflows to the database on hydrogen

"""

import sys
import pprint

import numpy as np

from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar
from pymatgen.io.vaspio.vasp_input import Potcar, Kpoints
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from fireworks.fw_config import LAUNCHPAD_LOC
from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket

from mpinterfaces import get_struct_from_mp
from mpinterfaces.firetasks import MPINTCalibrateTask
from mpinterfaces.firetasks import MPINTMeasurementTask
from mpinterfaces.firetasks import MPINTDatabaseTask
from mpinterfaces.interface import Ligand, Interface

def get_calibration_task(structure, phase="CalibrateBulk", \
                         slab_interface_params={'hkl':[1,0,0], 'ligand': None},\
                         turn_knobs={}, incar_params={}, other_params={}):
    """
    returns general calibration task for a structure
    
    Args:
        structure    : pymatgen structure to be calibrated (can be a bulk, ligand, slab
                       or interface)  
        phase        : calibration type, viz. CalibrateBulk, CalibrateMolecule,
                       CalibrateSlab, CalibrateInterface
        hkl          : in case of Slab and Interface miller indices of facet 
        turn_knobs   : specifies the parameters to be calibrated 
        incar_params : dictionary of additional incar parameters, refer defined 
                       incar_dict for defaults 
        other_params : other parameters for calibration, viz. job_dir, is_matrix, etc. 
                       described in the calibrate module
    """
    #structure definition 
    
    poscar = Poscar(structure)
    incar_dict = { 'SYSTEM': 'slab',
                   'ENCUT': 500, 
                   'ISIF': 2, 
                   'IBRION': 2, 
                   'ISMEAR': 1, 
                   'EDIFF': 1e-05, 
                   'NPAR': 4, 
                   'SIGMA': 0.1, 
                   'PREC': 'Accurate'
                 }
    if incar_params: 
        incar_dict.update(incar_params)
    incar = Incar.from_dict(incar_dict)
    kpoints = Kpoints.monkhorst_automatic(kpts=(8, 8, 1))
    que  = { 'nnodes':1,
             'nprocs':16,
             'walltime':'48:00:00',
             'job_bin': '/home/km468/Software/VASP/vaspsol_kappa.5.3.5/vasp'             
            }
    # calibration task: relax hkl
    calparams = {}
    calparams['calibrate'] = phase
    calparams['incar'] = incar.as_dict()
    calparams['poscar'] = poscar.as_dict()
    calparams['kpoints'] = kpoints.as_dict()
    calparams['que_params'] = que
    calparams['turn_knobs'] = turn_knobs
    if phase == 'CalibrateSlab':
         calparams['system'] = {'hkl':slab_interface_params['hkl'],
                                'ligand':slab_interface_params['ligand']
                               }
    elif phase == 'CalibrateInterface':
         calparams['system'] = {'hkl':hkl,
                                'ligand':structure.ligand.reduced_formula
                               }
    calparams['other_params'] = {
        'is_matrix':False,
        'from_ase':True,
        'Grid_type':'M'
        }
    if other_params:
        calparams['other_params'].update(other_params)
    return MPINTCalibrateTask(calparams)




def get_workflow(fireworks=None, name=None):
    """
    returns a workflow consisting of n fireworks
    with each firework dependent one on the other in sequence 
    """
    
    return Workflow(fireworks, links_dict = {fireworks[i-1].fw_id:fireworks[i].fw_id for \
            i,j in enumerate(fireworks)}, name=name )


def get_measurement_task(structure, sol_add={}): 
    """
    define a measurement task
    """ 

    # calibration task1: relax hkl
    #caltask = get_calibration_task(structure, hkl=hkl)
    # measurement task: static 
    msrparams = {}
    msrparams['measurement'] = 'MeasurementInterface'
    msrparams['que_params'] =  { 'nnodes':1,
                                 'nprocs':16,
                                 'walltime':'48:00:00',
                                 'job_bin': '/home/km468/Software/VASP/vaspsol_kappa.5.3.5/vasp'
                               }
    #msrparams['other_params'] = {'job_dir': name+'_static'}
    msrtask = MPINTMeasurementTask(msrparams)
    # measurement task: solvation
    solmsrparams = {}
    solmsrparams['measurement'] = 'MeasurementSolvation'
    solmsrparams['que_params'] =  msrparams['que_params']
    sol_params = { 'EB_K':[78.4],
                   'TAU':[],
                   'LAMBDA_D_K':[3.0],
                   'NELECT':list(np.linspace(-1,1,5))
                 }
    sol_params.update(sol_add)
    #solmsrparams['other_params'] = {'job_dir': name+'_sol',
    #                                'sol_params':sol_params
     #                              }
    return MPINTMeasurementTask(solmsrparams)
    #firework1 = [caltask]
    #workflow = [firework1, firework2, firework3]
    #firework2 is linked to firework1
    
if __name__=='__main__':
    """
    General workflow construction method: 
    1. construct firetasks
    2. collect firetasks 
    3. define fireworks 
    4. get workflow 
    for a ligand based set of workflows: 
    Read as 
    1. Workflow Name, within script structure pull
       1. Firework 1 (1., 2., Firetasks description)
       2. Firework 2 (1., 2., Firetasks description)

    1. CalBulk: first pull bulk from materialsproject or local POSCAR. 
       1. BulkConvergence: (1.cutoff convergence, kpoints convergence)
       2. BulkRelax: (1. relax the bulk)
    2. CalLigand: pull molecule from Smiles, construct Ligand
       1. LigandConvergence: cutoff convergence of Ligand
       2. LigandRelax: (1. relax the ligand) 
    3. CalSlab: Construct Slab from CONTCAR of BulkRelax
       1. SlabConvergence: Vacuum, Thickness convergence of Slab
    4. CalInterface: Construct Interface, Slab from CONTCARs of BulkRelax, LigandRelax,
       1. InterfaceConfig: (1.Configuration empirical relaxation) ###to be integrated
       2. InterfaceRelax: (1. Interface Relaxation of chosen POSCAR from firework 1)
       3. SlabRelax: (1. Slab Relaxation)
    5. Measurement
       1. SolLigand: (1. Solvation of Ligand)
    """




    ###################### FIRST WORKFLOW: BULK #######################################
    # pull structure, 
    # firework 1 firetask1: convergence 
    # firework 2 firetask1: relaxation

    workflows = []
    bulk_fireworks = [] 
    bulk_firetasks = []
    bulk = get_struct_from_mp('PbS', MAPI_KEY="2hnLK7D14uWUHoJs") #pull from matproj, use your own key 
#    primitive --> conventional cell
    sa = SpacegroupAnalyzer(bulk)
    structure_conventional = sa.get_conventional_standard_structure()
    bulk = structure_conventional.copy()
    bulk.sort()


    # create two calibration tasks: convergence of cutoff, kpoints and a relaxation
    bulk_firetasks.append(get_calibration_task(structure=bulk, \
                        turn_knobs= {'ENCUT': [400,500,600],
                                     'KPOINTS': [[x,x,x] for x in range(4,8)]}, \
                                     other_params={'job_dir': 'BulkConvergence'}))
    bulk_firetasks.append(get_calibration_task(structure=bulk, \
                                    incar_params= {'NSW':1000},\
                                    other_params={'job_dir': 'BulkRelax'}))
    
    # create fireworks 1. bulk_convergence, bulk_relax NOTE: update the 
    # NOTE: Update the second firework with incar param encut and kpoints based on 
    # desired convergence criteria 
    bulk_fireworks.append(Firework([bulk_firetasks[0]], name='bulk_convergence'))
    bulk_fireworks.append(Firework([bulk_firetasks[1]], name='bulk_relax'))

    # create workflows 
    workflows.append(get_workflow(bulk_fireworks, name='Bulk'))
    ######################################################################################




    ##################### SECOND WORKFLOW: LIGAND ########################################
    # pull structure of ligand 
    # firework 3 firetask1: cutoff convergence
    # firework 4 firetask1: relaxation  
    mols=[]
    lig_firetasks =[] 
    lig_fireworks = [] 
    molecule = Structure.from_file("POSCAR_TCPO")  #for testing purposes, can be constructed 
			                           #from smiles_to_poscar.py script
    mols.append(molecule)
    # create ligand from molecules
    ligand = Ligand(mols, cm_dist=[], angle={}, link={}, remove=[],
                 charge=0, spin_multiplicity=None,
                 validate_proximity=False)


    #create firetasks for molecule cutoff convergence and relaxation
    lig_firetasks.append(get_calibration_task(structure=ligand, \
                                              turn_knobs={'ENCUT':[400,500,600]}, \
                                              other_params={'job_dir': 'LigandConvergence'}))                        
    lig_firetasks.append(get_calibration_task(structure=ligand, \
                                    incar_params = {'NSW':100},\
                                    other_params={'job_dir': 'LigandRelax'}))

    #create firework for molecule cutoff convergence
    lig_fireworks.append(Firework(lig_firetasks[0], name='molecule_convergence'))


    workflows.append(get_workflow(lig_fireworks, name='Ligand'))
    #######################################################################################



    ################# THIRD WORKFLOW: SLAB Convergence #####################################
    # construct slab from relaxation firework CONTCAR of firework 2 
    # firework 4 firetask1: convergence of thickness and vacuum 
    slab_firetasks = [] 
    slab_fireworks = [] 
    #NOTE: This should be the CONTCAR file under the directory BulkRelax/CONTCAR
    #for testing purposes this is the same bulk structure pulled from matproj
    slab = Interface(bulk, hkl=[1,0,0], min_thick=10, min_vac=10,
                 supercell=[1,1,1], name=None, adsorb_on_species='Pb',
                 adatom_on_lig='O', ligand=None, displacement=3.0,
                 surface_coverage=0.01, scell_nmax=10, coverage_tol=0.25,
                 solvent="amine", start_from_slab=False, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False, primitive = True,
                 from_ase=True, x_shift= 0, y_shift= 0)

    slab_firetasks.append(get_calibration_task(structure=slab, \
                                              turn_knobs={'VACUUM':[400,500,600],
                                                          'THICKNESS':range(5,20)}))

    slab_fireworks.append(Firework(slab_firetasks[0], name='slab_convergence'))

    workflows.append(get_workflow(slab_fireworks, name='Slab'))

    #########################################################################################




    ################## FOURTH WORKFLOW: INTERFACE Configuration and SLAB Relax ##############
    # construct and configure interface from relaxed CONTCAR of 
    # firework 3 and firework 2 
    # at thickness and vac spacing of firework 4
    # firework 5 firetask1: Interface relaxation , firetask2: slab relaxation
    interface_firetasks = []
    interface_fireworks = []
    interface = Interface(bulk, hkl=[1,0,0], min_thick=10, min_vac=10,
                 supercell=[1,1,1], name=None, adsorb_on_species='Pb',
                 adatom_on_lig='O', ligand=ligand, displacement=3.0,
                 surface_coverage=0.01, scell_nmax=10, coverage_tol=0.25,
                 solvent="amine", start_from_slab=False, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False, primitive = True,
                 from_ase=True, x_shift= 0, y_shift= 0)
    
    interface.create_interface()  #with pre-configure as true, will generate different
                                  #ligand interface configurations
    interface.sort()
    slab = interface.slab
    slab.sort()
    
    # create interface firetasks: interface relaxation and slab relaxation 
    interface_firetasks.append(get_calibration_task(structure=interface, \
                                              incar_params = {'NSW':100}, \
                                              other_params = {'job_dir': 'InterfaceRelax'}))
    interface_firetasks.append(get_calibration_task(structure=slab, \
                                              incar_params = {'NSW':100}))

    interface_fireworks.append(Firework(interface_firetasks, name='interface_convergence'))

    workflows.append(get_workflow(interface_fireworks, name='Interface'))


    ##########################################################################################




    ######################### FIFTH WORKFLOW: Solvation ######################################
    # measurement fireworks: 3 firetasks for each slab, ligand , interface
    measurement_firetasks = []
    measurement_fireworks = []
    interface = interface
    slab =   slab #Structure.from_file("POSCAR")            #"from interface 4"
    ligand = ligand #Structure.from_file("POSCAR")            #"from interface 5"
    measurement_firetasks.append(get_measurement_task(structure=interface, \
                                              sol_add = {'EB_K':[1,4.6,47.6]}))
    measurement_firetasks.append(get_measurement_task(structure=slab, 
                                              sol_add = {'EB_K':[1,4.6,47.6]})) 
    measurement_firetasks.append(get_measurement_task(structure=ligand, 
                                              sol_add = {'EB_K': [1,4.6,47.6]})) 
    measurement_fireworks.append(Firework(measurement_firetasks, name='solvation'))
    workflows.append(get_workflow(measurement_fireworks, name='Solvation'))

    ############################Submit Workflows to Database#################################
    # connect to the fireworks database and add workflow to it
    # use your own account
    launchpad = LaunchPad.from_file(LAUNCHPAD_LOC)

    print('fireworks in the database before adding the workflow: \n',
        launchpad.get_fw_ids())

    for wf in workflows:
        launchpad.add_wf(wf, reassign_all=False)

    wfs = launchpad.get_wf_ids()
    pp = pprint.PrettyPrinter(indent=4)
    for i,fw_id in enumerate(wfs):
        print('Workflow: {}'.format(i+1))
        pp.pprint(launchpad.get_wf_summary_dict(fw_id))
        print('\n')
    print('Total number of workflows in the database: \n', len(wfs))
    print('Total number of fireworks in the database: \n', len(launchpad.get_fw_ids()))
        


