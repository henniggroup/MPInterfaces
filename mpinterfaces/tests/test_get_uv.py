""" tests for mpinterfaces.transformations """

from __future__ import unicode_literals

import unittest

from mpinterfaces.transformations import get_uv

__author__ = "Seve G. Monahan"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


class TestGetUV(unittest.TestCase):

    def test_uv_1_1(self):
        self.assertEqual([1, 1], get_uv([1, 1], [[1, 0], [0, 1]]))

    def test_uv_2_2(self):
        self.assertEqual([2, 2], get_uv([1, 1], [[2, 0], [0, 2]]))

    def test_uv_2_1(self):
        self.assertEqual([2, 1], get_uv([1, 1], [[2, 0], [0, 1]]))

    def test_uv_addition(self):
        self.assertEqual([5, 3], get_uv([2, 3], [[1, 1], [0, 1]]))

    def test_uv_subtraction(self):
        self.assertEqual([-1, 3], get_uv([2, 3], [[1, -1], [0, 1]]))

if __name__ == "__main__":
    unittest.main()
