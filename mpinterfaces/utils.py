from __future__ import division, unicode_literals, print_function

"""
Utility functions
"""

import sys
import math
import socket
import numpy as np

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from ase.lattice.surface import surface


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


def get_job_state(job_id):
    """
    return the job state given the job_id
    """
    hostname = socket.gethostname()
    state = None
    #hipergator,pbs
    if 'ufhpc' in hostname:
        output = sp.check_output(['qstat', '-i', job_id])
        state = output.rstrip('\n').split('\n')[-1].split()[-2]
    #stampede, slurm
    elif 'stampede' in hostname:
        output = sp.check_output(['squeue', '--job', job_id])
        state = output.rstrip('\n').split('\n')[-1].split()[-4]
    #no batch system
    else:
        state = 'XX'        
    return state
    
