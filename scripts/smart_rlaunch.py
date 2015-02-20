from __future__ import division, unicode_literals, print_function

import os
import mmap
from argparse import ArgumentParser

from fireworks.fw_config import LAUNCHPAD_LOC, FWORKER_LOC, CONFIG_FILE_DIR
from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.fworker import FWorker
from fireworks.core.rocket import Rocket

from fabric.api import settings, run

#MAILDIR = '/home/km468/.fireworks/mail/hipergator/'
MAILDIR = '/home/matk/Documents/mail/hipergator/'

def check_done(job_id):
    for fname in os.listdir(MAILDIR):
        #print fname, os.path.getmtime(MAILDIR+fname) - time.time()
        if present_id(job_id, MAILDIR+fname):
            print('id in mail {0} = {1}'.format(fname, job_id))            
            return True

def present_id(job_id, fname, search_string=b"Execution terminated"):
    with open(fname, 'r') as f:
        fm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        n = len(fm)
        if fm.rfind(search_string, 0, n) > 0:
            i = fm.rfind(b"PBS Job Id:", 0, n)
            j = fm.rfind(b"Job Name:", i, n)
            jid = fm[i:j].split()[-1]
            return job_id == jid
        else:
            return False
def launch():
    parser = ArgumentParser(description='smart rocket launcher')
    parser.add_argument('-f', '--fw_id', help='specific fw_id to run', default=None, type=int)
    args = parser.parse_args()
    job_ids = None
    lp = LaunchPad.from_file(LAUNCHPAD_LOC)
    m_fw = lp._get_a_fw_to_run(fw_id=args.fw_id, checkout=False)
    if m_fw is None:
        print('no firework with that id')
        return
    fw_spec = dict(m_fw.spec)
    #done = [False]*len(fw_spec['cal_objs'])
    done = []
    if fw_spec.get('cal_objs',None) is not None:
        for calparams in fw_spec['cal_objs']:
            if calparams.get('job_ids', None) is not None:
                job_ids = calparams.get('job_ids', None)
                print('job ids : 'job_ids)
                if job_ids is not None:
                    for jid in job_ids:
                        done.append(check_done(jid))
                else:
                    print('job_ids not set')
    else:
        print('this firework doesnt have any cal_objs')
        done.append(True)
    if done and all(done):
        with settings(host_string='km468@hipergator.rc.ufl.edu'):
            run("ssh dev1 rlaunch singleshot -f "+str(args.fw_id))
    else:
        print('Havent recieved execution termination confirmation from all the jobs in the firework')
        
if __name__ == '__main__':
    ##test
    #print(check_done('10995522.moab.ufhpc'))
    #print(check_done('10995523.moab.ufhpc'))
    #print(check_done('10995520.moab.ufhpc'))
    launch()

    
