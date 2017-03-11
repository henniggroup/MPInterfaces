""" tests for mpinterfaces.transformations """

from __future__ import unicode_literals

import unittest
from mpinterfaces.transformations import get_r_list
from io import StringIO
import sys

__author__ = "Seve G. Monahan"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


class TestGetRList(unittest.TestCase):

    def setUp(self):
        self.old_stdout = sys.stdout
        sys.stdout = self.mystdout = StringIO()

    def tearDown(self):
        sys.stdout = self.old_stdout

    # Expect rmax is the floor of the max divided by its
    # correponding area argument
    def test_half(self):
        self.assertEqual([[1, 1]], get_r_list(1, 0.98, 2))
        self.assertEqual('rmax1, rmax2: 2, 2\n\n', self.mystdout.getvalue())

    def test_half_decimal(self):
        self.assertEqual([[1, 1], [2, 2]], get_r_list(1.5, 1.5, 3))
        self.assertEqual('rmax1, rmax2: 2, 2\n\n', self.mystdout.getvalue())

    def test_half_error_margin(self):
        self.assertEqual([], get_r_list(1.5, 1.6, 3))
        self.assertEqual('rmax1, rmax2: 2, 1\n\n', self.mystdout.getvalue())

    def test_altered_error_margin(self):
        self.assertEqual([[1, 1]], get_r_list(1.5, 1.6, 3, 0.04))
        self.assertEqual('rmax1, rmax2: 2, 1\n\n', self.mystdout.getvalue())


if __name__ == "__main__":
    unittest.main()
