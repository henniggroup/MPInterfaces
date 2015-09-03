from __future__ import division, unicode_literals, print_function

import os
import sys
import time
from math import sqrt
from collections import OrderedDict

from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.inputs import Potcar, Kpoints
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from mpinterfaces import get_struct_from_mp
from mpinterfaces.interface import Interface
from mpinterfaces.calibrate import Calibrate, CalibrateSlab
from mpinterfaces.transformations import *
from mpinterfaces.utils import *


# default incar settings
incar_dict = dict(
    PREC = 'Accurate',
    ENCUT = 400,
    ISMEAR = 0,
    EDIFF = '1E-6',
    ISIF = 3,
    NSW = 500,
    IBRION = 2,
    NPAR = 4,
    LWAVE = '.FALSE.',
    LCHARG = '.FALSE.',
    GGA = 'BO',
    PARAM1 = 0.1833333333,
    PARAM2 = 0.2200000000,
    LUSE_VDW = '.TRUE.',
    AGGAC = 0.0000 )
incar_sub = Incar.from_dict(incar_dict)
incar_sub['ISMEAR'] = 1
incar_2d = Incar.from_dict(incar_dict)
# KPOINTS
kpoints_sub = Kpoints.monkhorst_automatic(kpts=(18, 18, 18))
kpoints_2d = Kpoints.monkhorst_automatic(kpts=(18, 18, 1))
# QUE
nprocs = 32
nnodes = 1
mem='1000'
walltime = '24:00:00'
bin_sub = '/home/km468/Software/VASP/vasp.5.3.5/vasp'
bin_2d = '/home/km468/Software/VASP/vasp.5.3.5/vasp_noz'
# STRUCTURES
substrates = [ 'Pt', 'Ag', 'Cu', 'Ni', 'Al' , 'Au', 'Pd', 'Ir']
mat2ds = ['POSCAR_graphene']


def run_cal(turn_knobs, qadapter, job_cmd, job_dir, name,
            incar=None, poscar=None, potcar=None, kpoints=None):
    """
    setup and run calibrate job
    """
    Calibrate.LOG_FILE = name+'.json'
    cal = Calibrate(incar, poscar, potcar, kpoints, 
                    turn_knobs=turn_knobs, qadapter=qadapter,
                    job_cmd = job_cmd, job_dir=job_dir)
    cal.setup()
    cal.run()

    
def step1():
    """
    get substrate bulk structures from materialsproject for Pt, Ag, Cu, Ni,
    Al, Au, Pd, Ir and do 3d relax
    
    get 2d structures from the provided poscars(just poscar_graphene)
    and relax in x and y only
    
    - POSCAR_graphene must be made avialable in the directory
    - creates required input files and submits the jobs to the que (9 jobs)
    - returns: step1_sub.json step1_2d.json 
    """
    #job directory for the runs
    job_dir_sub = 'step1_sub'
    job_dir_2d = 'step1_2d'
    # create list of all substrate poscars
    poscars_sub = []
    poscars_2d = []
    for sub in substrates:
        struct_sub = get_struct_from_mp(sub)
        sa_sub = SpacegroupAnalyzer(struct_sub)
        struct_sub = sa_sub.get_conventional_standard_structure()
        poscars_sub.append(Poscar(struct_sub))
    for td in mat2ds:
        poscars_2d.append(Poscar.from_file(td))
    # setup calibrate and run'em
    turn_knobs_sub = OrderedDict(
        [
            ('POSCAR', poscars_sub)
        ])
    turn_knobs_2d = OrderedDict(
        [
            ('POSCAR', poscars_2d)
        ])
    qadapter_sub, job_cmd_sub = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                            walltime=walltime,
                                            job_bin=bin_sub, mem=mem)
    qadapter_2d, job_cmd_2d = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                            walltime=walltime,
                                            job_bin=bin_2d, mem=mem)
    run_cal(turn_knobs_sub, qadapter_sub, job_cmd_sub, job_dir_sub,
            'step1_sub', incar=incar_sub, kpoints=kpoints_sub)
    run_cal(turn_knobs_2d, qadapter_2d, job_cmd_2d, job_dir_2d,
            'step1_2d', incar=incar_2d, kpoints=kpoints_2d)
    return ['step1_sub.json', 'step1_2d.json']


def step2():
    """
    read in relaxed bulk substrate, construct 111 slabs and relax only
    ionic positions, top 2 layers

    - uses info from step1_sub.json step1_2d.json
    - creates required input files and submits the jobs to the que	(9 jobs)
    - returns: step2.json    
    """
    hkl = [1,1,1]
    min_thick = 10
    min_vac = 18
    n_layers = 2
    #job directory for the runs
    name = 'step2'
    job_dir = 'step2'
    # incar & kpoints
    incar = Incar.from_dict(incar_dict)
    incar['ISMEAR'] = 1
    incar['IBRION'] = 2
    kpoints = Kpoints.monkhorst_automatic(kpts=(18, 18, 1))    
    # create list of all substrate poscars
    all_poscars = []
    # load in previous jobs
    all_jobs = Calibrate.jobs_from_file('step1_sub.json')    
    for job in all_jobs:
        job_dir = os.path.join(job.parent_job_dir, job.job_dir)
        contcar_file = os.path.join(job_dir, 'CONTCAR')
        relaxed_struct = Structure.from_file(contcar_file)
        substrate_slab = Interface(relaxed_struct,
                                hkl = hkl,
                                min_thick = min_thick,
                                min_vac = min_vac,
                                primitive = False, from_ase = True)
        sd_flags = CalibrateSlab.set_sd_flags(interface=substrate_slab,
                                            n_layers=n_layers)
        poscar = Poscar(substrate_slab, selective_dynamics=sd_flags)
        all_poscars.append(poscar)
    # setup calibrate and run'em
    turn_knobs = OrderedDict(
        [
            ('POSCAR', all_poscars)
        ])
    qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                      walltime=walltime,
                                      job_bin=bin_sub, mem=mem)
    run_cal(turn_knobs, qadapter, job_cmd, job_dir,
            name, incar=incar, kpoints=kpoints)
    return [name+'.json']


def step4():
    """
    put relaxed 2d materials in all possible ways on the relaxed slab and
   relax the interface structure, relax interface ionic positions

    merged step3 and 4
    - uses info from step2.json and step1_2d.json
    - creates required input files and submits the jobs to the que (around 40 jobs)   
    - returns: step4.json     
    """
    seperation = 3 # in angstroms
    nlayers_2d = 1
    nlayers_sub = 2
    hkl_sub = [1,1,1]
    hkl_2d = [0,0,1]
    #job directory for the runs
    name = 'step4'    
    job_dir = 'step4'
    # incar & kpoints
    incar = Incar.from_dict(incar_dict)
    incar['ISMEAR'] = 1
    incar['ISIF'] = 2
    kpoints = Kpoints.monkhorst_automatic(kpts=(18, 18, 1))
    # CSL settings for each substrate
    alignment_settings = { 'Pt': [120, 0.10, 1, 0.5],
                        'Ag': [120, 0.10, 1, 0.5],
                        'Al': [120, 0.10, 1, 0.5],
                        'Au': [120, 0.10, 1, 0.5],
                        'Pd': [120, 0.10, 1, 0.5],
                        'Ir': [120, 0.10, 1, 0.5],
                        'Cu': [50, 0.06, 1, 0.5],
                        'Ni': [50, 0.06, 1, 0.5] }
    # load in previous jobs
    relaxed_sub_jobs = Calibrate.jobs_from_file('step2.json')
    relaxed_2d_jobs = Calibrate.jobs_from_file('step1_2d.json')
    # create list of all substrate poscars
    all_poscars = []
    for jsub in relaxed_sub_jobs:
        job_dir_sub = os.path.join(jsub.parent_job_dir, jsub.job_dir)
        contcar_file = os.path.join(job_dir_sub, 'CONTCAR')
        relaxed_struct_sub = slab_from_file(hkl_sub, contcar_file)
        species_sub = ''.join([tos.symbol for tos in relaxed_struct_sub.types_of_specie])
        for j2d in relaxed_2d_jobs:    
            job_dir_2d = os.path.join(j2d.parent_job_dir, j2d.job_dir)
            contcar_file = os.path.join(job_dir_2d, 'CONTCAR')
            relaxed_struct_2d = slab_from_file(hkl_2d, contcar_file)
            species_2d = ''.join([tos.symbol for tos in relaxed_struct_2d.types_of_specie])
            print(species_sub, species_2d)
            substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(
                relaxed_struct_sub,
                relaxed_struct_2d,
                *alignment_settings[species_sub])
            sd_flags = CalibrateSlab.set_sd_flags(interface=mat2d_slab_aligned,
                                                n_layers=nlayers_2d)
            poscar = Poscar(mat2d_slab_aligned, selective_dynamics=sd_flags)
            poscar.comment = '_'.join([species_sub,species_2d,'2d'])
            all_poscars.append(poscar)
            sd_flags = CalibrateSlab.set_sd_flags(interface=substrate_slab_aligned,
                                                n_layers=nlayers_sub)
            poscar = Poscar(substrate_slab_aligned, selective_dynamics=sd_flags)
            poscar.comment = '_'.join([species_sub,species_2d,'sub'])
            all_poscars.append(poscar)        
            hetero_interfaces = generate_all_configs(mat2d_slab_aligned,
                                                    substrate_slab_aligned,
                                                    nlayers_2d,
                                                    nlayers_sub,
                                                    seperation )
            for i, iface in enumerate(hetero_interfaces):
                sd_flags = CalibrateSlab.set_sd_flags(
                    interface=iface,
                    n_layers=nlayers_2d+nlayers_sub,
                    top=True, bottom=False)
                poscar = Poscar(iface, selective_dynamics=sd_flags)
                poscar.comment = '_'.join([species_sub,species_2d,str(i)])
                all_poscars.append(poscar)
    # setup calibrate and run'em
    turn_knobs = OrderedDict(
        [
            ('POSCAR', all_poscars)
        ])
    qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                      walltime=walltime,
                                      job_bin=bin_sub, mem=mem)
    run_cal(turn_knobs, qadapter, job_cmd, job_dir,
            name, incar=incar, kpoints=kpoints)
    return [name+'.json']


def launch_daemon(steps):
    """
    run all th steps in daemon mode
    """
    for step in steps:
        while True:
            chkpt_files = globals()[step]()
            done = []            
            for cf in chkpt_files:
                time.sleep(3)                
                Calibrate.update_checkpoint(jfile=cf)
                all_jobs = Calibrate.jobs_from_file(cf)
                # test: done = [True]
                done = done + [(True if job.final_energy else False)
                               for job in all_jobs]
            if all(done):
                print('all jobs in {} done. Proceeding to the next one'.format(step))
                time.sleep(5)
                break
            # update and check every 5 mins
            time.sleep(300)
                    

if __name__ == '__main__':
    steps = ['step1', 'step2', 'step4']
    launch_daemon(steps)           
