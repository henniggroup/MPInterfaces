import os
import glob
from setuptools import setup, find_packages

SETUP_PTH = os.path.dirname(os.path.abspath(__file__))

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "mpinterfaces",
    version = "0.0.0",
    author = "Kiran Mathew, Joshua Gabriel",
    author_email = "km468@cornell.edu",
    description = ("vasp automation"),
    license = "GPL",
    url = "https://matk86@bitbucket.org/matk86/vasp_automation.git",
    packages=find_packages(),
    long_description=read('README'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: GPL License",
    ],
    scripts=glob.glob(os.path.join(SETUP_PTH, "scripts", "*"))
)
