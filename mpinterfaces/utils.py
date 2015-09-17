from __future__ import division, unicode_literals, print_function

"""
Utility functions
"""

import sys
import os
import math
import socket
import time
import subprocess as sp
import logging
import numpy as np

from monty.json import MSONable, MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.inputs import Potcar, Kpoints

from custodian.custodian import Custodian

from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from ase.lattice.surface import surface

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


def get_ase_slab(pmg_struct, hkl=(1,1,1), min_thick=10, min_vac=10):
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
    ase_slab.center(vacuum=min_vac/2, axis=2)
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
                Structure.from_sites(slab_input,to_unit_cell=True),
                shift=0,
                scale_factor=np.eye(3, dtype=np.int),
                site_properties=slab_input.site_properties)


def add_vacuum_padding(slab, vacuum, hkl=[0,0,1]):
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
    new_c = c / np.linalg.norm(c) * (thickness+vacuum)
    new_lattice = Lattice(np.array([a,b,new_c]))
    new_sites = []
    for site in slab:
        new_sites.append(PeriodicSite(site.species_and_occu,
                                      site.coords,
                                      new_lattice,
                                      properties=site.properties,
                                      coords_are_cartesian=True))
    new_struct = Structure.from_sites(new_sites)
    #center the slab
    avg_z = np.average([fcoord[2] for fcoord in new_struct.frac_coords])
    new_struct.translate_sites(list(range(len(new_struct))),
                               [0, 0, 0.5 - avg_z])
    return Slab(new_struct.lattice,
                new_struct.species_and_occu,
                new_struct.frac_coords,
                hkl,
                Structure.from_sites(new_struct,to_unit_cell=True),
                shift=0,
                scale_factor=np.eye(3, dtype=np.int),
                site_properties=new_struct.site_properties)


def get_run_cmmnd(nnodes=1, nprocs=16, walltime='24:00:00',
                  job_bin=None, mem='1000'):
    d = {}
    job_cmd = None
    hostname = socket.gethostname()
    #hipergator
    if 'ufhpc' in hostname:
        if job_bin is None:
            job_bin='/home/km468/Software/VASP/vasp.5.3.5/vasp'
        else:
            job_bin = job_bin
        d = {'type':'PBS',
             'params':
                 {
                'nnodes': str(nnodes),
                'ppnode': str(int(nprocs/nnodes)),
                'walltime': walltime,
                'job_name': 'vasp_job',
                'email': 'mpinterfaces@gmail.com',
                'notification_options': 'ae',
                'pre_rocket': '#PBS -l pmem='+str(mem)+'mb',
                'rocket_launch': 'mpirun '+job_bin
                }
             }
    #stampede
    elif 'stampede' in hostname:
        if job_bin is None:
            job_bin='/home1/01682/km468/Software/VASP/vasp.5.3.5/vasp'
        else:
            job_bin = job_bin
        d = {'type':'SLURM',
             'params':
                 {
                'nodes': str(nnodes),
                'ntasks': str(nprocs),
                'walltime': walltime,
                'queue':'normal',
                'account':'TG-DMR050028N',
                'job_name': 'vasp_job',
                'rocket_launch': 'ibrun '+job_bin
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
        job_cmd=['ls', '-lt']
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
    #hipergator,pbs
    if 'ufhpc' in hostname:
        try:
            output = sp.check_output(['qstat', '-i', job.job_id])
            state = output.rstrip('\n').split('\n')[-1].split()[-2]
        except:
            print('Job {} not in the que'.format(job.job_id))
            state = "00" 
        ofname = "FW_job.out"
    #stampede, slurm
    elif 'stampede' in hostname:
        try:
            output = sp.check_output(['squeue', '--job', job.job_id])
            state = output.rstrip('\n').split('\n')[-1].split()[-4]
        except:
            print('Job {} not in the que'.format(job.job_id))
            state = "00"             
        ofname = "vasp_job-"+str(job.job_id)+".out"
    #no batch system
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
        job_ids: list of job ids
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
    qadater = None
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
            logger.info('setting job {0} in {1} to rerun'.format(j['job_id'], job.job_dir))
            contcar_file = job.job_dir+os.sep+'CONTCAR'
            poscar_file = job.job_dir+os.sep+'POSCAR'
            if os.path.isfile(contcar_file) and len(open(contcar_file).readlines()) != 0 :
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

    
def launch_daemon(steps, interval, handlers=None):
        """
        run all the 'steps' in daemon mode
        checks job status every 'interval' seconds
        also runs all the error handlers
        """
        for step in steps:
            chkpt_files = step()
            while True:
                done = []            
                for cf in chkpt_files:
                    time.sleep(3)                
                    update_checkpoint(jfile=cf)
                    all_jobs = jobs_from_file(cf)
                    for j in all_jobs:
                        state, ofname = get_job_state(j)
                        if j.final_energy:
                            done = done + [True]
                        elif state == 'R':
                            logger.info('job {} running'.format(j.job_id))
                            done = done + [False]
                        elif state in ['C', 'CF', 'F', '00']:
                            logger.error('Job {0} in {1} cancelled or failed. State = {2}'.format(j.job_id, j.job_dir,state))
                            done = done + [False]
                            if handlers:
                                logger.info('Investigating ... ')
                                os.chdir(j.job_dir)
                                for h in handlers:
                                    if ofname:
                                        h.output_filename = ofname
                                        if h.check():
                                            logger.error('Detected vasp errors {}'.format(h.errors))
                                os.chdir(j.parent_job_dir)
                        else:
                            logger.info('Job {0} pending. State = {1}'.format(j.job_id,state))
                            done = done + [False]
                # test:
                #done = [True, True]                            
                if all(done):
                    logger.info('all jobs in {} done. Proceeding to the next one'.format(step.func_name))                                    
                    time.sleep(5)
                    break
                logger.info('all jobs in {0} NOT done. Next update in {1} seconds'.format(step.func_name,interval))
                time.sleep(interval)
        
        
def get_convergence_data(jfile):
    """
    returns data dict in the following format
    {'Al':
          {'ENCUT': [ [500,1.232], [600,0.8798] ], 
            'KPOINTS':[ [], [] ]
          },
     'W': ...
    }

    Note: only ENCUT and KPOINTS
    """
    cutoff_jobs = jobs_from_file(jfile)
    data = {}
    for j in cutoff_jobs:
        jdir = os.path.join(j.parent_job_dir, j.job_dir)
        poscar_file = os.path.join(jdir, 'POSCAR')
        struct_m = Structure.from_file(poscar_file)
        species = ''.join([tos.symbol for tos in struct_m.types_of_specie])
        # energy / atom
        if data.get(species):
            data[species]['ENCUT'].append([j.vis.incar['ENCUT'],j.final_energy/len(struct_m)])
            data[species]['KPOINTS'].append([j.vis.kpoints.kpts,j.final_energy/len(struct_m)])
        else:
            data[species] = {}
            data[species]['ENCUT'] = []
            data[species]['KPOINTS'] = []
    return data    


def get_opt_params(data, species, param='ENCUT'):
    """
    return optimum parameter
    default: 1 meV/atom
    """
    sorted_list = sorted(data[species][param], key=lambda x:x[1])
    sorted_array = np.array(sorted_list)
    # 1meV/atom difference
    consecutive_diff = np.abs(sorted_array[:-1,1] - sorted_array[1:,1] - 0.001)
    min_index = np.argmin(consecutive_diff)
    return sorted_list[min_index][0]
    

def partition_jobs(turn_knobs, max_jobs):
    """
    divide turn_knobs into smaller turn_knobs so that each one of
    them has smaller max_jobs jobs
    """
    params_len = [len(v) for k,v in turn_knobs.items()]
    n_total_jobs =  reduce(lambda x,y:x*y, params_len)
    partition_size = int(n_total_jobs/max_jobs)
    max_index = np.argmax(params_len)
    max_len = max(params_len)
    max_key = turn_knobs.items()[max_index][0]
    partition = range(0,max_len,max(1,int(max_len/partition_size)))
    partition_1 = partition[1:]+[max_len]
    print('{0} list of length {1} will be partitioned into {2} chunks'.format(max_key,max_len,len(partition)))
    turn_knobs_list = []
    name_list = []
    for i,j in zip(partition,partition_1):
        ordered_list = []
        for k,v in turn_knobs.items():
            if k == max_key:
                tk_item = (k,v[i:j])
            else:
                tk_item = (k,v)
            ordered_list.append(tk_item)
        turn_knobs_list.append(OrderedDict(ordered_list))
        name_list.append('_'.join([str(i),str(j)]))
    return turn_knobs_list, name_list

        
def create_step(turn_knobs, qadapter, job_cmd, job_dir,
                name, incar=None, kpoints=None):
    """
    step wrapper
    """
    def step():
        run_cal(turn_knobs, qadapter, job_cmd, job_dir,
                name, incar=incar, kpoints=kpoints)
        return [name+'.json']
    return step
