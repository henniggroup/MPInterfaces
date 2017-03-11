"""
test for calibrate.py main_recipe

runs main_recipe of calibrate

and compares the calibrate.json with the test output json file

similar to test on interface
"""

from subprocess import call
import unittest
import json

__author__ = "Joshua J. Gabriel"
__copyright__ = "Copyright 2017, Henniggroup"
__version__ = "1.6"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

# This test is simply calling a shell script, which calls a python main recipe
# (Originally a function used for ad-hoc manual testing) and verifies it behaved
# correctly. The reason for using a python file to call a shell script is so
# automatic python testing tools, like nose2, will automatially run it.


class StartupTest(unittest.TestCase):

    def test_main_calibrate(self):
        out = call("./test_calibrate_main.sh")
        self.assertEqual(out, 0)
        test_json = json.load(open('calibrate.json', 'r'))
        ans_json = json.load(open('../../test_files/calibrate.json', 'r'))
        self.assertDictEqual(test_json, ans_json)


if __name__ == '__main__':
    unittest.main()
