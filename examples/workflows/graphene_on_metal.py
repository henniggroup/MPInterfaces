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


# default incar settings, with vdw
incar_dict = dict(
    PREC = 'Accurate',
    ENCUT = 400,
    ISMEAR = 0,
    EDIFF = '1E-6',
    ISIF = 3,
    IBRION = 2,
    NSW = 500,
    NPAR = 4,
    LCHARG = '.FALSE.',
    GGA = 'BO',
    PARAM1 = 0.1833333333,
    PARAM2 = 0.2200000000,
    LUSE_VDW = '.TRUE.',
    AGGAC = 0.0000 )
# INCAR
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
    get substrate bulk structures from materialsproject for
    Pt, Ag, Cu, Ni, Al, Au, Pd, Ir and do 3d relaxation(ISIF=3)
    
    get 2d structures from the provided poscars(just poscar_graphene)
    and relax in x and y only(vasp_noz bin)
    
    - POSCAR_graphene must be made available in the directory
    - creates required input files and submits the jobs to the que
    - 8 + 1 jobs
    - returns: step1_sub.json step1_2d.json 
    """
    #job directory for the runs
    job_dir_sub = 'step1_sub'
    job_dir_2d = 'step1_2d'
    # create list of all substrate poscars
    poscars_sub = []
    poscars_2d = []
    # substrate structures
    for sub in substrates:
        struct_sub = get_struct_from_mp(sub)
        sa_sub = SpacegroupAnalyzer(struct_sub)
        struct_sub = sa_sub.get_conventional_standard_structure()
        poscars_sub.append(Poscar(struct_sub))
    # 2d structures
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
    # normal binary
    qadapter_sub, job_cmd_sub = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                            walltime=walltime,
                                            job_bin=bin_sub, mem=mem)
    # binary with z constraint
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
    read in the realxed bulk substrates and relaxed 2d,
    create substrate slab,
    get aligned substrates and 2d,
    relax the aligned structures seperatly(only ionic positions, ISIF=2)

    - input from step1_sub.json and step1_2d.json
    - 8(pairs) * 2 = 16 jobs
    - returns step2.json
    """
    nlayers_2d = 1
    nlayers_sub = 2
    hkl_sub = [1,1,1]
    min_thick = 10.0
    min_vac = 18.0
    hkl_2d = [0,0,1]
    #job directory for the runs
    job_dir_sub = 'step2_sub'
    job_dir_2d = 'step2_2d'
    # isif = 2
    incar_sub['ISIF'] = 2
    incar_2d['ISIF'] = 2
    # kpoints
    kpoints_sub = Kpoints.monkhorst_automatic(kpts=(18, 18, 1))
    kpoints_2d = Kpoints.monkhorst_automatic(kpts=(18, 18, 1))    
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
    relaxed_sub_jobs = Calibrate.jobs_from_file('step1_sub.json')
    relaxed_2d_jobs = Calibrate.jobs_from_file('step1_2d.json')
    poscars_sub = []
    poscars_2d = []
    # create list of all aligned substrate and 2d slabs
    for jsub in relaxed_sub_jobs:
        jdir = os.path.join(jsub.parent_job_dir, jsub.job_dir)
        contcar_file = os.path.join(jdir, 'CONTCAR')
        relaxed_struct_sub = Structure.from_file(contcar_file)
        # create slab
        slab_sub = Interface(relaxed_struct_sub,  hkl = hkl_sub,
                             min_thick = min_thick, min_vac = min_vac,
                             primitive = False, from_ase = True)
        species_sub = ''.join([tos.symbol for tos in slab_sub.types_of_specie])
        # loop over 2d
        for j2d in relaxed_2d_jobs:    
            jdir = os.path.join(j2d.parent_job_dir, j2d.job_dir)
            contcar_file = os.path.join(jdir, 'CONTCAR')
            slab_2d = slab_from_file(hkl_2d, contcar_file)
            species_2d = ''.join([tos.symbol for tos in slab_2d.types_of_specie])
            print(species_sub, species_2d)
            # align
            slab_sub_aligned, slab_2d_aligned = get_aligned_lattices(
                slab_sub,
                slab_2d,
                *alignment_settings[species_sub])
            # aligned sub poscar
            sd_flags = CalibrateSlab.set_sd_flags(interface=slab_sub_aligned,
                                                  n_layers=nlayers_sub)
            poscar = Poscar(slab_sub_aligned, selective_dynamics=sd_flags)
            poscar.comment = '_'.join([species_sub,species_2d,'sub'])            
            poscars_sub.append(poscar)
            # aligned 2d slab
            sd_flags = CalibrateSlab.set_sd_flags(interface=slab_2d_aligned,
                                                  n_layers=nlayers_2d)
            poscar = Poscar(slab_2d_aligned, selective_dynamics=sd_flags)
            poscar.comment = '_'.join([species_sub,species_2d,'2d'])
            poscars_2d.append(poscar)
    # setup calibrate and run'em
    turn_knobs_sub = OrderedDict(
        [
            ('POSCAR', poscars_sub)
        ])
    turn_knobs_2d = OrderedDict(
        [
            ('POSCAR', poscars_2d)
        ])
    qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
                                            walltime=walltime,
                                            job_bin=bin_sub, mem=mem)
    run_cal(turn_knobs_sub, qadapter, job_cmd, job_dir_sub,
            'step2_sub', incar=incar_sub, kpoints=kpoints_sub)
    run_cal(turn_knobs_2d, qadapter, job_cmd, job_dir_2d,
            'step2_2d', incar=incar_2d, kpoints=kpoints_2d)
    return ['step2_sub.json', 'step2_2d.json']


def step3():
    """
    put aligned & relaxed 2d materials in all possible ways on the
    aligned & relaxed slab,
    relax interface ionic positions(ISIF=2)

    - uses info from step2_sub.json and step2_2d.json
    - creates required input files and submits the jobs to the que
    - 8(pairs) * 2(atoms in graphene basis) = 16 jobs
    - returns: step3.json     
    """
    seperation = 3 # in angstroms
    nlayers_2d = 1
    nlayers_sub = 2
    hkl_sub = [1,1,1]
    hkl_2d = [0,0,1]
    #job directory for the runs
    name = 'step3'    
    job_dir = 'step3'
    # incar
    incar = Incar.from_dict(incar_dict)
    incar['ISMEAR'] = 1
    incar['ISIF'] = 2
    # kpoints
    kpoints = Kpoints.monkhorst_automatic(kpts=(18, 18, 1))
    # load in previous jobs
    relaxed_sub_jobs = Calibrate.jobs_from_file('step2_sub.json')
    relaxed_2d_jobs = Calibrate.jobs_from_file('step2_2d.json')
    # create list of all substrate poscars
    all_poscars = []
    # loop over aligned & relaxed substrates and 2d
    for jsub, j2d in zip(relaxed_sub_jobs,relaxed_2d_jobs):
        # substrate
        job_dir_sub = os.path.join(jsub.parent_job_dir, jsub.job_dir)
        contcar_file = os.path.join(job_dir_sub, 'CONTCAR')
        # read in as structure object
        substrate_slab_aligned = Structure.from_file(contcar_file)
        species_sub = ''.join([tos.symbol for tos in substrate_slab_aligned.types_of_specie])
        # 2d
        job_dir_2d = os.path.join(j2d.parent_job_dir, j2d.job_dir)
        contcar_file = os.path.join(job_dir_2d, 'CONTCAR')
        # read in as structure object        
        mat2d_slab_aligned = Structure.from_file(contcar_file)
        species_2d = ''.join([tos.symbol for tos in mat2d_slab_aligned.types_of_specie])
        print(species_sub, species_2d)
        # position the aligned materials in all possible ways
        hetero_interfaces = generate_all_configs(mat2d_slab_aligned,
                                                 substrate_slab_aligned,
                                                 nlayers_2d,
                                                 nlayers_sub,
                                                 seperation )
        # loop over all hetero-interfaces
        for i, iface in enumerate(hetero_interfaces):
            sd_flags = CalibrateSlab.set_sd_flags(interface=iface,
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


def launch_daemon(steps, interval):
    """
    run all the steps in daemon mode
    checks job status every 'interval' seconds
    """
    for step in steps:
        chkpt_files = globals()[step]()        
        while True:
            done = []            
            for cf in chkpt_files:
                time.sleep(3)                
                Calibrate.update_checkpoint(jfile=cf)
                all_jobs = Calibrate.jobs_from_file(cf)
                # test:
                #done = [True, True]
                done = done + [(True if job.final_energy else False)
                               for job in all_jobs]
            if all(done):
                print('all jobs in {} done. Proceeding to the next one'.format(step))
                time.sleep(5)
                break
            # update and check every 'interval' seconds
            print('all jobs in {0} NOT done. Next update in {1} seconds'.format(step,interval))
            time.sleep(interval)
                    

if __name__ == '__main__':
    # name of functions to be run.
    # functions will be run in the order given in the list
    steps = ['step1', 'step2', 'step3']
    # update interval
    interval = 300
    launch_daemon(steps,interval)
