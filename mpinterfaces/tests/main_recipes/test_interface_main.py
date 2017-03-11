from subprocess import call

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
    assert(call("./test_interface_main.sh") == 0)

if __name__ == '__main__':
    test_main_interface()
