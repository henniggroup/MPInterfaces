from subprocess import call
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.structure import Structure
from glob import glob
import os 

# This test is simply calling a shell script, which calls a python main recipe
# (Originally a function used for ad-hoc manual testing) and verifies it behaved
# correctly. The reason for using a python file to call a shell script is so
# automatic python testing tools, like nose2, will automatially run it.

__author__ = "Seve G. Monahan"
__copyright__ = "Copyright 2017, Henniggroup"
__version__ = "1.6"
__maintainer__ = "Joshua J. Gabriel"
__email__ = "joshgabriel92@gmail.com"
__status__ = "Production"
__date__ = "March 3, 2017"

def test_main_interface():
    # check screen output
    assert(call("./test_interface_main.sh") == 0)
    # check intermediate and output files
    for f in glob('*.vasp'):
        #print (f)
        try:
           s1 = Structure.from_file('{}'.format(f))
           s2 = Structure.from_file('../../test_files/{}'.format(f))
           assert(StructureMatcher().fit(s1,s2))
        except:
           print (f)

    

if __name__ == '__main__':
    test_main_interface()
    os.system('rm POSCAR_diacetate_boxed.vasp POSCAR_interface.vasp POSCAR_slab.vasp lead_acetate.xyz')
