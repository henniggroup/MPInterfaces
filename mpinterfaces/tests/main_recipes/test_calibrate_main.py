"""
test for calibrate.py main_recipe

runs main_recipe of calibrate

and compares the calibrate.json with the test output json file

similar to test on interface
"""
from subprocess import call
import pandas as pd
# This test is simply calling a shell script, which calls a python main recipe
# (Originally a function used for ad-hoc manual testing) and verifies it behaved
# correctly. The reason for using a python file to call a shell script is so
# automatic python testing tools, like nose2, will automatially run it.

def test_main_calibrate():
    assert(call("./test_calibrate_main.sh") == 0)
    test_json  = pd.read_json('calibrate.json')
    assert(test_json.equals(pd.read_json('../../test_files/calibrate.json')))

if __name__ == '__main__':
    test_main_calibrate()
