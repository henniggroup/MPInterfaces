# basic test for instrument

from __future__ import unicode_literals

import unittest

import os

import json

from mpinterfaces.instrument import *
from pymatgen.io.vasp.inputs import Incar,Kpoints,Poscar,Potcar

__author__ = "Joshua J. Gabriel"
__copyright__ = "Copyright 2017, Henniggroup"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

TEST_STEP1 = os.path.join(os.path.dirname(__file__), "..", "test_files/Wflow_Step1")
TEST_STEP2 = os.path.join(os.path.dirname(__file__), "..", "test_files/Wflow_Step2")

class TestInstrument(unittest.TestCase):

    def test_write_inputset(self):
        name = 'Test'
        incar= Incar.from_file(TEST_STEP1+os.sep+'INCAR')
        kpoints = Kpoints.from_file(TEST_STEP1+os.sep+'KPOINTS')
        poscar = Poscar.from_file(TEST_STEP1+os.sep+'POSCAR')
        potcar = TEST_STEP1+os.sep+'DUMMY_POTSPEC'
        #potcar = #Potcar.from_dict({'@class': 'Potcar', 'functional': 'PBE',\
                  #                 'symbols': ['Al'], '@module': 'pymatgen.io.vasp.inputs'})
        reuse_path = [TEST_STEP1 + os.sep + 'COPY_FILE']
        print (reuse_path)
        mvis = MPINTVaspInputSet(name,incar,poscar,potcar,kpoints,reuse_path=reuse_path,test=True)
        mvis.write_input(job_dir=TEST_STEP2)
        self.assertCountEqual(os.listdir(TEST_STEP2), ['INCAR','KPOINTS','POSCAR','COPY_FILE'])
        cleanup = [os.remove(TEST_STEP2+os.sep+f) for f in os.listdir(TEST_STEP2)]

if __name__ == '__main__':
    TI = TestInstrument()
    TI.test_write_inputset()
