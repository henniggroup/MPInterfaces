#!/usr/bin/env python
"""
fetch a single firework from the database and run it on the machine where this script
is invoked
"""

from fireworks import LaunchPad
from fireworks.core.rocket_launcher import launch_rocket #rapidfire

if __name__ = '__main__':
    launchpad = LaunchPad(host='localhost', port=27017, name='fireworks',
                        username="km468", password="km468")

    print 'fireworks in the database before adding the workflow: \n', launchpad.get_fw_ids()
    launch_rocket(launchpad)#, fworker=None, fw_id=None, strm_lvl='INFO')
    print 'fireworks in the database: \n', launchpad.get_fw_ids()
    #rapidfire(launchpad)
