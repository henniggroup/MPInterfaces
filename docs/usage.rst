Launchpad queries
==================

Connecting to the Database
---------------------------

The mongo database for job submission('fireworks') is set up on the
machine 'hydrogen'.
Please use your own account to connect to the databse
Contact me(km468@cornell.edu) to create a database account

Since hydrogen is part of the  uf network any machine onthe uf
network(inlcuding hipergator) should be able to access hydrogen
at 10.5.46.101
        
To connect to the database from outside uf network( for example from
stampede), tunnel port number 27017 from that machine to port 27017
on hydrogen via ssh:

   ssh -N -f -L 27017:10.5.46.101:27017 username@hipergator.rc.ufl.edu

if port 27017 on the machine that you are running is not available,
use another port number for tunneling. example:-

   ssh -N -f -L 27030:10.5.46.101:27017 username@hipergator.rc.ufl.edu

mind: if the tunneled port is changed, the port number in the
launchpad initialization should also be changed



Workflows
----------

- general info:
      to submit workflow to the database:

         python workflow_script.py
 
      Fireworks package has some nice utility scripts for launching
      fireworks and checking job status. If the fireworks package is
      installed then those scripts are already in your PATH. Some
      examples are given below

      initialize database connection(writes a yaml file with the 
      database settings in the directory where it is called).
      Tip: dont have to do this again if you create a ~/.fireworks
      folder and the generated yaml file there:

         lpad init

	 CAUTION: Be careful when using the following command as it will 
	 erase all workflows from the database:

         lpad reset

	 launch a single firework:

         rlaunch singleshot

	 get all fireworks info:

         lpad get_fws

- query by name:
      example:
	lpad get_wflows -n "Ag_100"

- query by state
      lpad get_wflows -s STATE
      
      or
      
      lpad get_wflows -q '{"state":STATE}' --sort created_on

      where STATE = "READY", "WAITING", "RUNNING", "FIZZLED" or "DEFUSED"

- delete workflows:
      example:
          lpad delete_wflows -n Ag_100
    

Fireworks
----------

- query by state
      lpad get_fws -s STATE
      
      or
      
      lpad get_fws -q '{"state":STATE}' --sort created_on

      where STATE = "READY", "WAITING", "RUNNING", "FIZZLED", "DEFUSED"

- query fireworks by name or id:
      example:
           lpad get_fws -n "solvation"
	   
	   lpad get_fws -i 102 -d all

- re-run firework with id, fw_id. Same as marking the firework as ready
     lpad rerun_fws -i fw_id

- re-run a finished or fizzled firework with updated specs:
       example: update the sol_params of the first task of the firework with id 102
		lpad rerun_fws -i 102
		
  		lpad update_fws -i 102 -u '{"_tasks.0.other_params.sol_params.NELECT":[-1,-0.5,0,0.5,1]}'
