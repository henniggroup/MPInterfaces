# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
Utility functions
"""

from six.moves import range
from six.moves import zip
from functools import reduce

import sys
import os
import math
import socket
import time
import subprocess as sp
import logging
import numpy as np
from collections import OrderedDict

from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Poscar

from custodian.custodian import Custodian

from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from ase.lattice.surface import surface

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


def get_ase_slab(pmg_struct, hkl=(1, 1, 1), min_thick=10, min_vac=10):
    """
    takes in the intial structure as pymatgen Structure object
    uses ase to generate the slab
    returns pymatgen Slab object
    
    Args:
        pmg_struct: pymatgen structure object
        hkl: hkl index of surface of slab to be created
        min_thick: minimum thickness of slab in Angstroms
        min_vac: minimum vacuum spacing 
    """
    ase_atoms = AseAtomsAdaptor().get_atoms(pmg_struct)
    pmg_slab_gen = SlabGenerator(pmg_struct, hkl, min_thick, min_vac)
    h = pmg_slab_gen._proj_height
    nlayers = int(math.ceil(pmg_slab_gen.min_slab_size / h))
    ase_slab = surface(ase_atoms, hkl, nlayers)
    ase_slab.center(vacuum=min_vac / 2, axis=2)
    pmg_slab_structure = AseAtomsAdaptor().get_structure(ase_slab)
    return Slab(lattice=pmg_slab_structure.lattice,
                species=pmg_slab_structure.species_and_occu,
                coords=pmg_slab_structure.frac_coords,
                site_properties=pmg_slab_structure.site_properties,
                miller_index=hkl, oriented_unit_cell=pmg_slab_structure,
                shift=0., scale_factor=None, energy=None)


def slab_from_file(hkl, filename):
    """
    reads in structure from the file and returns slab object.
    useful for reading in 2d/substrate structures from file.
    Args:
         hkl: miller index of the slab in the input file.
         filename: structure file in any format 
                   supported by pymatgen
    Returns:
         Slab object
    """
    slab_input = Structure.from_file(filename)
    return Slab(slab_input.lattice,
                slab_input.species_and_occu,
                slab_input.frac_coords,
                hkl,
                Structure.from_sites(slab_input, to_unit_cell=True),
                shift=0,
                scale_factor=np.eye(3, dtype=np.int),
                site_properties=slab_input.site_properties)


def add_vacuum_padding(slab, vacuum, hkl=[0, 0, 1]):
    """
    add vacuum spacing to the given structure
    Args:
        slab: sructure/slab object to be padded 
        vacuum: in angstroms
        hkl: miller index
    Returns:
         Structure object
    """
    min_z = np.min([fcoord[2] for fcoord in slab.frac_coords])
    slab.translate_sites(list(range(len(slab))), [0, 0, -min_z])
    a, b, c = slab.lattice.matrix
    z = [coord[2] for coord in slab.cart_coords]
    zmax = np.max(z)
    zmin = np.min(z)
    thickness = zmax - zmin
    new_c = c / np.linalg.norm(c) * (thickness + vacuum)
    new_lattice = Lattice(np.array([a, b, new_c]))
    new_sites = []
    for site in slab:
        new_sites.append(PeriodicSite(site.species_and_occu,
                                      site.coords,
                                      new_lattice,
                                      properties=site.properties,
                                      coords_are_cartesian=True))
    new_struct = Structure.from_sites(new_sites)
    # center the slab
    avg_z = np.average([fcoord[2] for fcoord in new_struct.frac_coords])
    new_struct.translate_sites(list(range(len(new_struct))),
                               [0, 0, 0.5 - avg_z])
    return Slab(new_struct.lattice,
                new_struct.species_and_occu,
                new_struct.frac_coords,
                hkl,
                Structure.from_sites(new_struct, to_unit_cell=True),
                shift=0,
                scale_factor=np.eye(3, dtype=np.int),
                site_properties=new_struct.site_properties)


<<<<<<< HEAD
def get_run_cmmnd(nnodes=1, nprocs=16, walltime='24:00:00',
                  job_bin=None, mem='1000'):
    d = {}
    job_cmd = None
    hostname = socket.gethostname()
    # hipergator
    if 'ufhpc' in hostname:
        if job_bin is None:
            job_bin = '/home/km468/Software/VASP/vasp.5.3.5/vasp'
        else:
            job_bin = job_bin
        d = {'type': 'PBS',
             'params':
                 {
                     'nnodes': str(nnodes),
                     'ppnode': str(int(nprocs / nnodes)),
                     'walltime': walltime,
                     'job_name': 'vasp_job',
                     'email': 'mpinterfaces@gmail.com',
                     'notification_options': 'ae',
                     'pre_rocket': '#PBS -l pmem=' + str(mem) + 'mb',
                     'rocket_launch': 'mpirun ' + job_bin
=======
def get_run_cmmnd(nnodes=1, ntasks=16, walltime='24:00:00',
                  job_bin=None, mem='1000', job_name=None):
    """
    depends on the supercomputing faciltiy being used.
    set a sample submit script in the fireworks directory which is
    installed in your virtual environment as for example:
    /my_venv/lib/python2.7/site-packages/FireWorks-1.2.5-py2.7.egg/
    fireworks/user_objects/queue_adapters/
    the keys to the dictionary d are the defaults on ufhpc's
    hipergator2 supercomputing facility

    """
    d = {}
    job_cmd = None
    hostname = socket.gethostname()
    # hipergator: currently hipergator2
    if 'ufhpc' in hostname:
        if job_bin is None:
            job_bin = '/home/mashton/vasp.5.4.1/bin/vasp'
        else:
            job_bin = job_bin
        d = {'type': 'SLURM',
             'params':
                 {
                     'nodes': str(nnodes),
                     'ntasks': str(int(ntasks)),
                     'walltime': walltime,
                     'job_name': job_name,
                     'email': 'mpinterfaces@gmail.com',
                     'notification_options': 'ae',
                     'pre_rocket': 'module load intel/2016.0.109 openmpi',
                     'rocket_launch': 'mpiexec ' + job_bin
>>>>>>> bb_real/master
                 }
             }
    # stampede
    elif 'stampede' in hostname:
        if job_bin is None:
            job_bin = '/home1/01682/km468/Software/VASP/vasp.5.3.5/vasp'
        else:
            job_bin = job_bin
        d = {'type': 'SLURM',
             'params':
                 {
                     'nodes': str(nnodes),
                     'ntasks': str(nprocs),
                     'walltime': walltime,
                     'queue': 'normal',
                     'account': 'TG-DMR050028N',
                     'job_name': 'vasp_job',
                     'rocket_launch': 'ibrun ' + job_bin
                 }
             }
    # running henniggroup machines
    elif hostname in ['hydrogen', 'helium',
                      'lithium', 'beryllium',
                      'carbon']:
        job_cmd = ['nohup', '/opt/openmpi_intel/bin/mpirun',
                   '-n', str(nprocs),
                   job_bin]
    # test
    else:
        job_cmd = ['ls', '-lt']
    if d:
        return (CommonAdapter(d['type'], **d['params']), job_cmd)
    else:
        return (None, job_cmd)


def get_job_state(job):
    """
    Args:
        job: job

    Returns:
           the job state and the job output file name
    """
    hostname = socket.gethostname()
    state = None
    ofname = None
    # hipergator,pbs
    if 'ufhpc' in hostname:
        try:
            output = sp.check_output(['qstat', '-i', job.job_id])
            state = output.rstrip('\n').split('\n')[-1].split()[-2]
        except:
            logger.info('Job {} not in the que'.format(job.job_id))
            state = "00"
        ofname = "FW_job.out"
    # stampede, slurm
    elif 'stampede' in hostname:
        try:
            output = sp.check_output(['squeue', '--job', job.job_id])
            state = output.rstrip('\n').split('\n')[-1].split()[-4]
        except:
            logger.info('Job {} not in the que.'.format(job.job_id))
            logger.info(
                'This could mean either the batchsystem crashed(highly unlikely) or the job completed a long time ago')
            state = "00"
        ofname = "vasp_job-" + str(job.job_id) + ".out"
    # no batch system
    else:
        state = 'XX'
    return state, ofname


def update_checkpoint(job_ids=None, jfile=None, **kwargs):
    """
    rerun the jobs with job ids in the job_ids list. The jobs are
    read from the json checkpoint file, jfile. 
    If no job_ids are given then the checkpoint file will 
    be updated with corresponding final energy
    Args:
        job_ids: list of job ids to update or q resolve
        jfile: check point file
    """
    cal_log = loadfn(jfile, cls=MontyDecoder)
    cal_log_new = []
    all_jobs = []
    run_jobs = []
    handlers = []
    final_energy = None
    incar = None
    kpoints = None
    qadapter = None
    # if updating the specs of the job
    for k, v in kwargs.items():
        if k == 'incar':
            incar = v
        if k == 'kpoints':
            kpoints = v
        if k == 'que':
            qadapter = v
    for j in cal_log:
        job = j["job"]
        job.job_id = j['job_id']
        all_jobs.append(job)
        if job_ids and (j['job_id'] in job_ids or job.job_dir in job_ids):
            logger.info('setting job {0} in {1} to rerun'.format(j['job_id'],
                                                                 job.job_dir))
            contcar_file = job.job_dir + os.sep + 'CONTCAR'
            poscar_file = job.job_dir + os.sep + 'POSCAR'
            if os.path.isfile(contcar_file) and len(
                    open(contcar_file).readlines()) != 0:
                logger.info('setting poscar file from {}'
                            .format(contcar_file))
                job.vis.poscar = Poscar.from_file(contcar_file)
            else:
                logger.info('setting poscar file from {}'
                            .format(poscar_file))
                job.vis.poscar = Poscar.from_file(poscar_file)
            if incar:
                logger.info('incar overridden')
                job.vis.incar = incar
            if kpoints:
                logger.info('kpoints overridden')
                job.vis.kpoints = kpoints
            if qadapter:
                logger.info('qadapter overridden')
                job.vis.qadapter = qadapter
            run_jobs.append(job)
    if run_jobs:
        c = Custodian(handlers, run_jobs, max_errors=5)
        c.run()
    for j in all_jobs:
        final_energy = j.get_final_energy()
        cal_log_new.append({"job": j.as_dict(),
                            'job_id': j.job_id,
                            "corrections": [],
                            'final_energy': final_energy})
    dumpfn(cal_log_new, jfile, cls=MontyEncoder,
           indent=4)


def jobs_from_file(filename='calibrate.json'):
    """
    read in json file of format caibrate.json(the default logfile
    created when jobs are run through calibrate) and return the
    list of job objects.
    Args:
        filename: checkpoint file name
    Returns:
           list of all jobs
    """
    caljobs = loadfn(filename, cls=MontyDecoder)
    all_jobs = []
    for j in caljobs:
        job = j["job"]
        job.job_id = j['job_id']
        job.final_energy = j['final_energy']
        all_jobs.append(job)
    return all_jobs


def launch_daemon(steps, interval, handlers=None, ld_logger=None):
    """
    run all the 'steps' in daemon mode
    checks job status every 'interval' seconds
    also runs all the error handlers
    """
    if ld_logger:
        global logger
        logger = ld_logger
    chkpt_files_prev = None
    for step in steps:
        chkpt_files = step(checkpoint_files=chkpt_files_prev)
        chkpt_files_prev = chkpt_files
        if not chkpt_files:
            return None
        while True:
            done = []
            reruns = []
            for cf in chkpt_files:
                time.sleep(3)
                update_checkpoint(job_ids=reruns, jfile=cf)
                all_jobs = jobs_from_file(cf)
                for j in all_jobs:
                    state, ofname = get_job_state(j)
                    if j.final_energy:
                        done = done + [True]
                    elif state == 'R':
                        logger.info('job {} running'.format(j.job_id))
                        done = done + [False]
                    elif state in ['C', 'CF', 'F', '00']:
                        logger.error(
                            'Job {0} in {1} cancelled or failed. State = {2}'. \
                                format(j.job_id, j.job_dir, state))
                        done = done + [False]
                        if handlers:
                            logger.info('Investigating ... ')
                            os.chdir(j.job_dir)
                            if ofname:
                                if os.path.exists(ofname):
                                    for h in handlers:
                                        h.output_filename = ofname
                                        if h.check():
                                            logger.error(
                                                'Detected vasp errors {}'.format(
                                                    h.errors))
                                            # TODO: correct the error and mark the job for rerun
                                            # all error handling must done using proper errorhandlers
                                            # h.correct()
                                            # reruns.append(j.job_id)
                                else:
                                    logger.error(
                                        'stdout redirect file not generated, job {} will be rerun'.format(
                                            j.job_id))
                                    reruns.append(j.job_id)
                            os.chdir(j.parent_job_dir)
                    else:
                        logger.info(
                            'Job {0} pending. State = {1}'.format(j.job_id,
                                                                  state))
                        done = done + [False]
            if all(done):
                logger.info(
                    'all jobs in {} done. Proceeding to the next one'.format(
                        step.__name__))
                time.sleep(5)
                break
            logger.info(
                'all jobs in {0} NOT done. Next update in {1} seconds'.format(
                    step.__name__, interval))
            time.sleep(interval)


def get_convergence_data(jfile, params=['ENCUT', 'KPOINTS']):
    """
    returns data dict in the following format
    {'Al':
          {'ENCUT': [ [500,1.232], [600,0.8798] ], 
            'KPOINTS':[ [], [] ]
          },
     'W': ...
    }

    Note: processes only INCAR parmaters and KPOINTS
    """
    cutoff_jobs = jobs_from_file(jfile)
    data = {}
    for j in cutoff_jobs:
        jdir = os.path.join(j.parent_job_dir, j.job_dir)
        poscar_file = os.path.join(jdir, 'POSCAR')
        struct_m = Structure.from_file(poscar_file)
        species = ''.join([tos.symbol for tos in struct_m.types_of_specie])
        if data.get(species):
            for p in params:
                if j.vis.incar.get(p):
                    data[species][p].append([j.vis.incar[p],
                                             j.final_energy / len(struct_m)])
                elif p == 'KPOINTS':
                    data[species]['KPOINTS'].append([j.vis.kpoints.kpts,
                                                     j.final_energy / len(
                                                         struct_m)])
                else:
                    logger.warn(
                        'dont know how to parse the parameter {}'.format(p))
        else:
            data[species] = {}
            for p in params:
                data[species][p] = []
                data[species][p] = []
    return data


def get_opt_params(data, species, param='ENCUT', ev_per_atom=0.001):
    """
    return optimum parameter
    default: 1 meV/atom
    """
    sorted_list = sorted(data[species][param], key=lambda x: x[1])
    sorted_array = np.array(sorted_list)
    consecutive_diff = np.abs(
        sorted_array[:-1, 1] - sorted_array[1:, 1] - ev_per_atom)
    min_index = np.argmin(consecutive_diff)
    return sorted_list[min_index][0]


# PLEASE DONT CHANGE THINGS WITHOUT UPDATING SCRIPTS/MODULES THAT DEPEND
# ON IT
# get_convergence_data and get_opt_params moved to *_custom
def get_convergence_data_custom(jfile, params=['ENCUT', 'KPOINTS']):
    """
    returns data dict in the following format
    {'Al':
          {'ENCUT': [ [500,1.232], [600,0.8798] ], 
            'KPOINTS':[ [], [] ]
          },
     'W': ...
    }

    Note: processes only INCAR parmaters and KPOINTS
    Sufficient tagging of the data assumed from species,Poscar
    comment line and potcar functional
    """
    cutoff_jobs = jobs_from_file(jfile)
    data = {}
    for j in cutoff_jobs:
        jdir = os.path.join(j.parent_job_dir, j.job_dir)
        poscar_file = os.path.join(jdir, 'POSCAR')
        struct_m = Structure.from_file(poscar_file)

        species = ''.join([tos.symbol for tos in struct_m.types_of_specie])
        tag = '_'.join([species, Poscar.from_file(poscar_file).comment,
                        j.vis.potcar.functional])
        if data.get(tag):
            for p in params:
                if j.vis.incar.get(p):
                    data[tag][p].append([j.vis.incar[p],
                                         j.final_energy / len(struct_m),
                                         j.vis.potcar, j.vis.poscar])
                    #                print(j.vis.potcar.functional,j.vis.poscar)
                elif p == 'KPOINTS':
                    data[tag]['KPOINTS'].append([j.vis.kpoints.kpts,
                                                 j.final_energy / len(
                                                     struct_m), j.vis.potcar,
                                                 j.vis.poscar])
                else:
                    logger.warn(
                        'dont know how to parse the parameter {}'.format(p))
        else:
            data[tag] = {}
            for p in params:
                data[tag][p] = []
                data[tag][p] = []
    return data


def get_opt_params_custom(data, tag, param='ENCUT', ev_per_atom=1.0):
    """
    Args:
        data:  dictionary of convergence data
        tag:   key to dictionary of convergence dara
        param: parameter to be optimized
        ev_per_atom: minimizing criterion in eV per unit 
    
    Returns
        [list] optimum parameter set consisting of tag, potcar object,
        poscar object, list of convergence data energies sorted according to 
        param
    
    default criterion: 1 meV/atom
    """
    sorted_list = sorted(data[tag][param], key=lambda x: x[0])
    # sorted array data
    t = np.array(sorted_list)[:, 1]
    # print(sorted_array[:-1,1], sorted_array[1:,1], ev_per_atom)
    consecutive_diff = [float(j) - float(i) - ev_per_atom for i, j in
                        zip(t[:-1], t[1:])]
    # print("Consecutive_diff",consecutive_diff)
    min_index = np.argmin(consecutive_diff)
    # return the tag,potcar object, poscar object, incar setting and convergence data for plotting that is optimum
    return [tag, data[tag][param][min_index][2],
            data[tag][param][min_index][3], sorted_list[min_index][0], t]


def partition_jobs(turn_knobs, max_jobs):
    """
    divide turn_knobs into smaller turn_knobs so that each one of
    them has smaller max_jobs jobs
    """
    params_len = [len(v) for k, v in turn_knobs.items()]
    n_total_jobs = reduce(lambda x, y: x * y, params_len)
    partition_size = int(n_total_jobs / max_jobs)
    max_index = np.argmax(params_len)
    max_len = max(params_len)
    max_key = list(turn_knobs.items())[max_index][0]
    partition = range(0, max_len, max(1, int(max_len / partition_size)))
    partition_1 = partition[1:] + [max_len]
    logger.info(
        '{0} list of length {1} will be partitioned into {2} chunks'.format(
            max_key, max_len, len(partition)))
    turn_knobs_list = []
    name_list = []
    for i, j in zip(partition, partition_1):
        ordered_list = []
        for k, v in turn_knobs.items():
            if k == max_key:
                tk_item = (k, v[i:j])
            else:
                tk_item = (k, v)
            ordered_list.append(tk_item)
        turn_knobs_list.append(OrderedDict(ordered_list))
        name_list.append('_'.join([str(i), str(j)]))
    return turn_knobs_list, name_list


def get_logger(log_file_name):
    """
    writes out logging file. 
    Very useful project logging, recommended for use 
    to monitor the start and completion of steps in the workflow
    Arg:
        log_file_name: name of the log file, log_file_name.log
    """
    loggr = logging.getLogger(log_file_name)
    loggr.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    fh = logging.FileHandler(log_file_name + '.log', mode='a')
    fh.setFormatter(formatter)
    loggr.addHandler(fh)
    return loggr


def set_sd_flags(poscar_input=None, n_layers=2, top=True, bottom=True,
                 poscar_output='POSCAR2'):
    """
    set the relaxation flags for top and bottom layers of interface.        
    The upper and lower bounds of the z coordinate are determined
    based on the slab. 
    Args:
         poscar_input: input poscar file name
         n_layers: number of layers to be relaxed
         top: whether n_layers from top are be relaxed
         bottom: whether n_layers from bottom are be relaxed
         poscar_output: output poscar file name
    Returns:
         None
         writes the modified poscar file
    """
    poscar1 = Poscar.from_file(poscar_input)
    sd_flags = np.zeros_like(poscar1.structure.frac_coords)
    z_coords = poscar1.structure.frac_coords[:, 2]
    z_lower_bound, z_upper_bound = None, None
    if bottom:
        z_lower_bound = np.unique(z_coords)[n_layers - 1]
        sd_flags[np.where(z_coords <= z_lower_bound)] = np.ones((1, 3))
    if top:
        z_upper_bound = np.unique(z_coords)[-n_layers]
        sd_flags[np.where(z_coords >= z_upper_bound)] = np.ones((1, 3))
    poscar2 = Poscar(poscar1.structure, selective_dynamics=sd_flags.tolist())
    poscar2.write_file(filename=poscar_output)
