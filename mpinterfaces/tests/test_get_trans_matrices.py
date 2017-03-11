""" tests for mpinterfaces.transformations """

from __future__ import unicode_literals

import unittest
from mpinterfaces.transformations import get_trans_matrices

__author__ = "Seve G. Monahan"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

# We're going to turn the generator into a list,
# because it really should be.


def list_trans(number):
    """ Turn the get_trans_matrices of number into a list, with each pair
        of numbers being a single entry on that list. The is preceding a
        refactor where get_trans_matrices will output the output of this
        function """
    generated_list = list(get_trans_matrices(number))

    generated_list_2 = []

    for i in generated_list:
        for j in i:
            generated_list_2.append(j)

    return generated_list_2


class TestGetTransMatrices(unittest.TestCase):
    """ Test get_trans_matrices through list_trans """

    def test_0(self):
        self.assertEqual([], list_trans(0))

    def test_1(self):
        self.assertEqual([[[1, 0], [0, 1]]], list_trans(1))

    def test_2(self):
        self.assertEqual([[[1, 0], [0, 2]],
                          [[1, 1], [0, 2]],
                          [[2, 0], [0, 1]]], list_trans(2))

    def test_3(self):
        self.assertEqual([[[1, 0], [0, 3]],
                          [[1, 1], [0, 3]],
                          [[1, 2], [0, 3]],
                          [[3, 0], [0, 1]]], list_trans(3))

    def test_9(self):
        self.assertEqual([[[1, 0], [0, 9]],
                          [[1, 1], [0, 9]],
                          [[1, 2], [0, 9]],
                          [[1, 3], [0, 9]],
                          [[1, 4], [0, 9]],
                          [[1, 5], [0, 9]],
                          [[1, 6], [0, 9]],
                          [[1, 7], [0, 9]],
                          [[1, 8], [0, 9]],
                          [[3, 0], [0, 3]],
                          [[3, 1], [0, 3]],
                          [[3, 2], [0, 3]],
                          [[9, 0], [0, 1]]], list_trans(9))

    def test_10(self):
        self.assertEqual([[[1, 0], [0, 10]],
                          [[1, 1], [0, 10]],
                          [[1, 2], [0, 10]],
                          [[1, 3], [0, 10]],
                          [[1, 4], [0, 10]],
                          [[1, 5], [0, 10]],
                          [[1, 6], [0, 10]],
                          [[1, 7], [0, 10]],
                          [[1, 8], [0, 10]],
                          [[1, 9], [0, 10]],
                          [[2, 0], [0, 5]],
                          [[2, 1], [0, 5]],
                          [[2, 2], [0, 5]],
                          [[2, 3], [0, 5]],
                          [[2, 4], [0, 5]],
                          [[5, 0], [0, 2]],
                          [[5, 1], [0, 2]],
                          [[10, 0], [0, 1]]], list_trans(10))

    def test_12(self):
        self.assertEqual([[[1, 0], [0, 12]],
                          [[1, 1], [0, 12]],
                          [[1, 2], [0, 12]],
                          [[1, 3], [0, 12]],
                          [[1, 4], [0, 12]],
                          [[1, 5], [0, 12]],
                          [[1, 6], [0, 12]],
                          [[1, 7], [0, 12]],
                          [[1, 8], [0, 12]],
                          [[1, 9], [0, 12]],
                          [[1, 10], [0, 12]],
                          [[1, 11], [0, 12]],
                          [[2, 0], [0, 6]],
                          [[2, 1], [0, 6]],
                          [[2, 2], [0, 6]],
                          [[2, 3], [0, 6]],
                          [[2, 4], [0, 6]],
                          [[2, 5], [0, 6]],
                          [[3, 0], [0, 4]],
                          [[3, 1], [0, 4]],
                          [[3, 2], [0, 4]],
                          [[3, 3], [0, 4]],
                          [[4, 0], [0, 3]],
                          [[4, 1], [0, 3]],
                          [[4, 2], [0, 3]],
                          [[6, 0], [0, 2]],
                          [[6, 1], [0, 2]],
                          [[12, 0], [0, 1]]], list_trans(12))

if __name__ == "__main__":
    unittest.main()
