""" tests for mpinterfaces """

from __future__ import unicode_literals

import unittest

from mpinterfaces.default_logger \
import get_default_logger

from io import StringIO

__author__ = "Seve G. Monahan"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"


class TestGetDefaultLogger(unittest.TestCase):
    def setUp(self):
        self.mystdout = StringIO()
        self.logger = get_default_logger("Example_Module_name", self.mystdout)

    def test_log(self):
        self.logger.info('Test')
        self.assertEqual('INFO:Example_Module_name:Test\n',\
                self.mystdout.getvalue())

if __name__ == "__main__":
    unittest.main()
