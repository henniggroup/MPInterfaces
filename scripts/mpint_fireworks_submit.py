"""
Create a workflow and submit it to the 'fireworks' database on hydrogen

Note 1:
      use own account to run the script. contact me(km468@cornell.edu) to create a database account
Note 2:
        forward port number 27017 from your local machine to port 27017 on hydrogen:
        ssh -N -f -L 27017:10.1.255.101:27017 username@hermes.mse.ufl.edu
"""

import numpy as np
from pymatgen import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Potcar, Kpoints
from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket #rapidfire
from mpinterfaces.firetasks import MPINTCalibrateTask


# set up the LaunchPad
launchpad = LaunchPad(host='localhost', port=27017, name='fireworks',
                       username="km468", password="km468")

# create the individual FireWorks and Workflow
system = 'Pt bulk'
atoms = ['Pt']
    
a0 = 3.965
lvec = [ [0.5, 0.0, 0.5], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5] ]
lvec = np.array(lvec) * a0
lattice = Lattice(lvec)#.from_parameters(3.866, 3.866, 3.866, 60, 60, 60)
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

#set the paramters for the calibrate task
calparams = {}
#required params: incar, poscar, list of encut and kpoints and the calibrate object
calparams['incar'] = incar.as_dict()
calparams['poscar'] = poscar.as_dict()
#range specification for encut and kpoints
calparams['encut_list'] = ['400', '800', '100']
calparams['kpoint_list'] = ['[7,7,7]', '[11,11,11]' ]
#type of calibration to be done: basically the name of calibrate calss to
#be used. available options: CalibrateMolecule, CalibrateSlab, CalibrateBulk
calparams['calibrate'] = 'CalibrateMolecule'
#optional param: job_dir is the name of the directory within which
#the encut and kpoints jobs will be run
calparams['cal_construct_params'] = {'job_dir':'Molecule_1'}

#create a calibrate task 
caltask1 = MPINTCalibrateTask(calparams)

calparams['calibrate'] = 'CalibrateSlab'
#optional param: job_dir is the name of the directory within which
#the encut and kpoints jobs will be run
calparams['cal_construct_params'] = {'job_dir':'Slab_1'}

#create a calibrate task 
caltask2 = MPINTCalibrateTask(calparams)



fw1 = Firework([caltask1, caltask2], name="calibrate")
wf = Workflow([fw1], name="mpint workflow")
#fw2 = Firework(measuretask, name="measurement", parents=[fw1])
#fw3 = Firework(pptask, name="post_process", parents=[fw1, fw2])
#wf = Workflow([fw1, fw2, fw3], name="mpint workflow")

# add workflow to the database
print 'fireworks in the database before adding the workflow: \n', launchpad.get_fw_ids()
launchpad.add_wf(wf)
print 'fireworks in the database: \n', launchpad.get_fw_ids()

