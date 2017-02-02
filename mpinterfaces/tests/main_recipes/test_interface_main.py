from subprocess import call

# This test is simply calling a shell script, which calls a python main recipe
# (Originally a function used for ad-hoc manual testing) and verifies it behaved
# correctly. The reason for using a python file to call a shell script is so
# automatic python testing tools, like nose2, will automatially run it.

def test_main_interface():
    assert(call("./test_interface_main.sh") == 0)

if __name__ == '__main__':
    test_main_interface()
