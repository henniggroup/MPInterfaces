from fireworks import LaunchPad
from fireworks.core.rocket_launcher import launch_rocket #rapidfire


# set up the LaunchPad
launchpad = LaunchPad(host='localhost', port=27017, name='fireworks',
                       username="km468", password="km468")
#reset launchpad
# Create a new FireWorks database. This will overwrite the existing FireWorks database!
#To safeguard against accidentally erasing an existing database, a password must
#be entered.
#
#launchpad.reset('', require_password=False)

launch_rocket(launchpad)#, fworker=None, fw_id=None, strm_lvl='INFO')
print launchpad.get_fw_ids()
#rapidfire(launchpad)
