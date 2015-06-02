import os
import glob
from setuptools import setup, find_packages

MPINT_DIR = os.path.dirname(os.path.abspath(__file__))

setup(
    name = "mpinterfaces",
    version = "0.1",
    install_requires=["numpy", "scipy",
                      "pymatgen", "custodian",
                      "fireworks", "pymatgen-db"],
    extras_require={"plotting": ["matplotlib>=1.1"],
                    "ase": ["ase>=3.3"],
                    "babel": ["openbabel", "pybel"],
                    "remote_launch": ["fabric"],
                    "documentation": ["sphinx"]
                   },    
    author = "Kiran Mathew, Joshua Gabriel",
    author_email = "km468@cornell.edu",
    description = ("High throughput analysis of interfaces using VASP and Materials Project tools"),
    license = "GPL",
    url = "https://matk86@bitbucket.org/matk86/vasp_automation.git",
    packages=find_packages(),
    long_description=open(os.path.join(os.path.dirname(__file__), 'README.rst')).read(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Development Status :: Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GPL License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
#    scripts=glob.glob(os.path.join(MPINT_DIR, "scripts", "*"))
)
