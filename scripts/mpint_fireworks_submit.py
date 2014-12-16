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
#reset launchpad
# Create a new FireWorks database. This will overwrite the existing FireWorks database!
#To safeguard against accidentally erasing an existing database, a password must
#be entered.
#
#launchpad.reset('', require_password=False)

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
calparams['incar'] = incar.as_dict()
calparams['poscar'] = poscar.as_dict()
#range specification for encut and kpoints
calparams['encut_list'] = ['400', '800', '100']
calparams['kpoint_list'] = ['[7,7,7]', '[11,11,11]' ]

caltask = MPINTCalibrateTask(calparams)
fw1 = Firework(caltask, name="calibrate")
wf = Workflow([fw1], name="mpint workflow")
#fw2 = Firework(measuretask, name="measurement", parents=[fw1])
#fw3 = Firework(pptask, name="post_process", parents=[fw1, fw2])
#wf = Workflow([fw1, fw2, fw3], name="mpint workflow")

# add workflow
launchpad.add_wf(wf)
#launch_rocket(launchpad)#, fworker=None, fw_id=None, strm_lvl='INFO')
print launchpad.get_fw_ids()
#rapidfire(launchpad)
