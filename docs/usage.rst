Overview
========

The MPInterfaces package contains python modules that enable the creation of 
FireWorks workflows for ligand-slab interface and slab-slab interface
calculations with VASP.

The package is structured as python modules (listed in modules section 
of docs and contained in mpinterfaces folder), examples notebooks 
(ipython notebooks) and workflows scripts to be 
used with a locally set up FireWorks database. 

The examples folder contains scripts for interface creation and utility post processing and 
ipython notebook examples and workflow submission scripts. 

Apart from these, there are utility scripts for the plotting of band structures, 
density of states and equations of state. 

The example folder "notebooks" also provides example calculations of Pourbaix diagrams, 
and Phase diagrams based on the method implemented in pymatgen/examples

Capabilities
============

MPInterfaces extends object oriented definitions of materialsproject pymatgen to 
interfaces and features:  

Interface Creation 
------------------
extended object oriented definitions to 
the pymatgen Slab called "Interface" and an 
adsorbate extension of pymatgen Molecule 
called "Ligand" and scripts to create Interfaces. 

Example scripts: 
    - create_interface.py : creates an example ligand + slab interface illustrating
    - ligand_interface.py : creates a list of ligand + slab interfaces with pre-configuration
    - hetero_interface.py : creates a slab – slab heterostructure interface 
    - nanoparticle.py     : creates a nanoparticle based on Wulff Construction


Offline vasp project management
--------------------------------

The script, **mpint**, in the scripts folder serves as a management
tool for vasp projects, starting from encut, kpoint or other parameter
optimization till the solvation calculations. Just define all types of
calculations with their corresponding specifications needed for the
project in a yaml file and run or rerun calculaitons as required.

it takes 3 arguments: input yaml file, type of calculation and the
run mode

example:

   ```
   mpint -i naf.yaml -t bulk_calibrate run
   ```
   
   this will read in the specifications for 'bulk_calibrate' job
   from the input yaml file, naf.yaml(in examples folder), and
   runs the job i.e submits to the PBS que.
   
Everytime jobs are submitted or its status queried, information
such as job ids, job folders etc are written to the calibrate.json 
checkpoint file. This makes it easier to identify job ids and their
corresponding job folders. Also provides all job info including 
converged final energies.in one place.

update the calibrate.json checkpoint file with the final energies(if converged):

       mpint update

update calibrate.json and rerun jobs with the given ids:

       mpint update 14692739.moab.ufhpc 14692740.moab.ufhpc

update calibrate.json and rerun jobs in the given folders:

       mpint update all_poscars/POS/CBr2_294_C1Br2 all_poscars/POS/CoBr2_294_Co1Br2

update calibrate.json and rerun jobs in the given folders with incar or kpoints or que parameters overridden by the parmaters in the override_input.yaml file:
       mpint -i override_input.yaml update all_poscars/POS/CBr2_294_C1Br2 all_poscars/POS/CoBr2_294_Co1Br2

High throughput VASP workflows with FireWorks  
---------------------------------------------

High throughput VASP workflows can be handled by creating your own FireWorks database server. 
The FireWorks package from the materialsproject is installed at the time of installation. 
We have defined FireTasks customized to Interface studies as follows: 

Calibrate, Measurement and DatabaseSubmit. 
    -  Calibrate Task: 
       Defines all tasks that involve the calibration of  
       VASP input parameters. For example, energy cutoff convergence, 
       kpoints convergence,vacuum convergence, thickness convergence in case of slabs 
       Calibrate module defined in calibrate.py provides customized 
       Bulk, Molecule, Slab and Interface objects for the Calibrate Task

    -  Measurement Task: 
       Tasks for VASP jobs dependent on a previous calculation. Currently implemented 
       tasks include Static calculation following a Calibrate Task and Solvation calculation following a Calibrate Task 
       Measurement tasks are linked to calibrate tasks and are capable of copying 
       the necessary files from the previous calculation to restart the calculation.
       for example, CONTCAR from a Relaxation for a Static calculation, CONTCAR and WAVECAR for a Solvation calculation

    -  Database Task: 
       Submits data to a MongoDB database similar to VasptoDBtaskDrone 
       of pymatgen extended to use Interface specific keys of miller index of slab 
       and ligand chemical formula in .json format

Refer examples/Workflows for example workflow submission scripts. 

Example of FireWorks Database Queries
=====================================

After installation and setup of your database according to 
FireWorks documentation, create your own my_launchpad.yaml using

lpad init

Set the path to this file to LAUNCHPAD_LOC in 
the fw_config.py found in your installation directory.  

TIP: create the my_launchpad.yaml in the directory
     ~/.fireworks

Make sure to configure your my_launchpad.yaml file to your local 
database settings (host name, database name, etc.)


Workflow submission and launching
---------------------------------

Refer examples folder for simple workflow scripts

      -> general submission of workflow to the database:

      python workflow_script.py
 
      Fireworks package has some nice utility scripts for launching
      fireworks and checking job status. If the fireworks package is
      installed then those scripts are already in your PATH. Some
      examples are given below:

      -> reset database to erase all workflows

	 CAUTION: Be careful when using the following command as it will 
	 erase all workflows from the database:

         lpad reset

      -> launch a single firework:

         rlaunch singleshot -f [fw_id]

      -> get all fireworks info:

         lpad get_fws

      -> get all workflows info:

         lpad get_wflows

      -> reset database to erase all workflows

         CAUTION: Be careful when using the following command as it will
         erase all workflows from the database:

         lpad reset



Workflow query examples
------------------------

      -> query by name:
           
           example:
	   lpad get_wflows -n "Ag_100"

      -> query by state
           
           lpad get_wflows -s STATE
      
           or
      
           lpad get_wflows -q '{"state":STATE}' --sort created_on

           where STATE = "READY", "WAITING", "RUNNING", "FIZZLED" or "DEFUSED"

      -> delete workflows:
           
           example:
           lpad delete_wflows -n Ag_100
    

Query for individual Fireworks
------------------------------

      -> query by state
           lpad get_fws -s STATE
      
      or
      
           lpad get_fws -q '{"state":STATE}' --sort created_on

           where STATE = "READY", "WAITING", "RUNNING", "FIZZLED", "DEFUSED"

      -> query fireworks by name or id:
           
           example:
           lpad get_fws -n "solvation"
	   
	   lpad get_fws -i 102 -d all

      -> re-run firework with id, fw_id. Same as marking the firework as ready
           lpad rerun_fws -i fw_id

      -> re-run a finished or fizzled firework with updated specs:
           example: update the sol_params of the first task of the firework with id 102
	   lpad rerun_fws -i 102
		
           lpad update_fws -i 102 -u '{"_tasks.0.other_params.sol_params.NELECT":[-1,-0.5,0,0.5,1]}'


Connecting to our local Database
--------------------------------

The mongo database for job submission('fireworks') is set up on the
machine 'hydrogen'.
Please use your own account to connect to the database
Contact me (km468@cornell.edu) to create a database account

