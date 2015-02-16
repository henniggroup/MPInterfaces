Workflows
==========
- query workflows by name:
      example:
	lpad get_wflows -n Ag_100

- query by state
      lpad get_wflows -s STATE
      
      or
      
      lpad get_wflows -q '{"state":STATE}' --sort created_on

      where STATE = "READY", "WAITING", "RUNNING", "FIZZLED", "DEFUSED"

- delete workflows:
      example:
          lpad delete_wflows -n Ag_100
    

Fireworks
==========
- query by state
      lpad get_fws -s STATE
      
      or
      
      lpad get_fws -q '{"state":STATE}' --sort created_on

      where STATE = "READY", "WAITING", "RUNNING", "FIZZLED", "DEFUSED"

- query fireworks by name or id:
      example:
           lpad get_fws -n "solvation"
	   
	   lpad get_fws -i 102 -d all

- re-run the firework id=fw_id(same as marking the firework with as ready)
     lpad rerun_fws -i fw_id

- re-run a finished or fizzled firework with updated specs:
       example: update the sol_params of the first task of the firework with id 102
		lpad rerun_fws -i 102
		
  		lpad update_fws -i 102 -u '{"_tasks.0.other_params.sol_params.NELECT":[-1,-0.5,0,0.5,1]}'
