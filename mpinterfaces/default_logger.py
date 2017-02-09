"""
Creates a default logger for outputing text.
"""
from __future__ import division, print_function, unicode_literals, \
        absolute_import

import os
import logging

def get_default_logger(name):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(formatter)
    logger.addHandler(sh)


