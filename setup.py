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
    description = ("Automate vasp calculations for interfaces"),
    license = "GPL",
    url = "https://matk86@bitbucket.org/matk86/vasp_automation.git",
    packages=find_packages(),
    long_description=read('README'),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GPL License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    scripts=glob.glob(os.path.join(SETUP_PTH, "scripts", "*"))
)
