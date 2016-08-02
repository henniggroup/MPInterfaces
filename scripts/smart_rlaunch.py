from __future__ import division, unicode_literals, print_function, \
    absolute_import

from six.moves import range

import mmap
import os
import time
from argparse import ArgumentParser

from fabric.api import settings, run
from fireworks import LaunchPad
from fireworks.fw_config import LAUNCHPAD_LOC

MAILDIR = os.path.join(os.path.expanduser('~'), ".fireworks/mail/hipergator/")


def check_done(job_id):
    for fname in os.listdir(MAILDIR):
        # print fname, os.path.getmtime(MAILDIR+fname) - time.time()
        if present_id(job_id, MAILDIR + fname):
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
    descr = "Smart rocket launcher." \
            "Uses the execution termination emails send by the batch system to " \
            "launch fireworks that depend on other fireworks i.e fireworks that " \
            "have 'cal_objs' in their spec." \
            "Takes in fw_ids of fireworks to be launched. " \
            "Range specification of fw_ids is also supported." \
            "All the jobs are launched to hipergator remotely using fabric." \
            "Note: Ensure that the submit scripts have the appropriate email settings and " \
            "that the mail daemon is fetching emails to the folder polled by this script." \
            "Note: All rocket launches take place in the home directory. If " \
            "_launch_dir spec is set for the firework the job files will be written to " \
            "that folder and jobs to the batch system will be done from there."
    parser = ArgumentParser(description=descr)
    parser.add_argument('-f', '--fw_ids', help='one or more of fw_ids to run',
                        default=None, type=int, nargs='+')
    parser.add_argument('-r', '--r_fw_ids',
                        help='start and end fw_ids of the range of fw_ids to run',
                        default=None, type=int, nargs=2)
    args = parser.parse_args()
    fw_ids = args.fw_ids
    if args.r_fw_ids is not None:
        fw_ids += list(range(args.r_fw_ids[0], args.r_fw_ids[1]))
    job_ids = None
    lp = LaunchPad.from_file(LAUNCHPAD_LOC)
    print('Firework ids: ', fw_ids)
    if fw_ids is None:
        print('No fw ids given')
        return
    for fid in fw_ids:
        m_fw = lp._get_a_fw_to_run(fw_id=fid, checkout=False)
        if m_fw is None:
            print('No firework with that id')
            return
        fw_spec = dict(m_fw.spec)
        done = []
        if fw_spec.get('cal_objs', None) is not None:
            for calparams in fw_spec['cal_objs']:
                if calparams.get('job_ids', None) is not None:
                    job_ids = calparams.get('job_ids', None)
                    print(fid, ' depends on jobs with ids : ', job_ids)
                    if job_ids is not None:
                        for jid in job_ids:
                            done.append(check_done(jid))
                    else:
                        print('job_ids not set')
        else:
            print('This firework doesnt depend on any other fireworks.')
            done.append(True)
        if done and all(done):
            print('Launching ', fid, ' ...')
            with settings(host_string='username@hipergator.rc.ufl.edu'):
                run("ssh dev1 rlaunch singleshot -f " + str(fid))
        else:
            print(
                "Haven't recieved execution termination confirmation for the jobs in the firework from hipergator resource manager")
        time.sleep(3)
    return


if __name__ == '__main__':
    launch()
    ##test
    # print(check_done('10995522.moab.ufhpc'))
    # print(check_done('10995523.moab.ufhpc'))
    # print(check_done('10995520.moab.ufhpc'))
