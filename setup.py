import os
import glob
from setuptools import setup, find_packages

MPINT_DIR = os.path.dirname(os.path.abspath(__file__))

setup(
    name = "mpinterfaces",
    version = "1.1.1",
    install_requires=["numpy>=1.8", "scipy>=0.10",
                      "pymatgen>=3.1.0", "fireworks>=1.07",
                      "custodian>=0.8.2", "pymatgen-db>=0.4.6",
                      "ase>=3.9.0"],
    extras_require={"plotting": ["matplotlib>=1.4.2"],
                    "babel": ["openbabel", "pybel"],
                    "remote_launch": ["fabric"],
                    "documentation": ["sphinx>=1.3.1", "sphinx-rtd-theme>=0.1.8"]
                   },    
    author = "Kiran Mathew, Joshua Gabriel, Richard G. Hennig",
    author_email = "km468@cornell.edu",
    description = ("High throughput analysis of interfaces using VASP and Materials Project tools"),
    license = "GPL",
    url = "https://matk86@bitbucket.org/matk86/vasp_automation.git",
    packages=find_packages(),
    long_description=open(os.path.join(os.path.dirname(__file__), 'README.rst')).read(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    scripts=glob.glob(os.path.join(MPINT_DIR, "scripts", "*"))
)
