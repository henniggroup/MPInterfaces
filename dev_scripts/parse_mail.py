import mmap
import os
import time

# MAILDIR = '/home/km468/.fireworks/mail/hipergator/'
MAILDIR = '/home/matk/Documents/mail/hipergator/'


def get_id(fname):
    with open(fname, 'r') as f:
        fm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        n = len(fm)
        i = fm.rfind(b"PBS Job Id:", 0, n)
        j = fm.rfind(b"Job Name:", i, n)
        return fm[i:j].split()[-1]


if __name__ == '__main__':
    for fname in os.listdir(MAILDIR):
        print fname, os.path.getmtime(MAILDIR + fname) - time.time()
        print 'id', get_id(MAILDIR + fname)
