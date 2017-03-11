"""
Creates a default logger for outputing text.
"""
from __future__ import division, print_function, unicode_literals, \
        absolute_import

import logging
import sys


def get_default_logger(name, output_stream=sys.stdout):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
    sh = logging.StreamHandler(stream=output_stream)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger


