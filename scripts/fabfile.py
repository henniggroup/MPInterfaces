from __future__ import with_statement

"""
remote job launch using the fabric package
usage:
    example: launch firework with id 10 on lithium
        fab launch:10,lithium

Note:
    the job will be launched from the home directory, unless
    an explicit launch directory is given using the
    firework spec '_launch_dir'(see workflow_examples)
"""

from fabric.api import settings, run, env, cd
from fabric.contrib.console import confirm

def launch(fw_id, machine='hipergator'):
    cmd = "rlaunch singleshot -f "+str(fw_id)    
    if machine == 'hipergator':
        with settings(host_string='km468@hipergator.rc.ufl.edu'):
            run("ssh dev1 "+cmd)
    elif machine == 'hydrogen':
        with settings(host_string='km468@hermes.mse.ufl.edu'):
            run("ssh hydrogen "+cmd)
    elif machine == 'stampede':
        pass
    else:
        with settings(host_string='km468@hermes.mse.ufl.edu'):
            run("ssh hydrogen ssh "+machine+" "+cmd)
            
if __name__ == '__main__':
    launch()
