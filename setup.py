import os
import glob
from setuptools import setup, find_packages

MPINT_DIR = os.path.dirname(os.path.abspath(__file__))

setup(
    name = "mpinterfaces",
    version = "1.1.2",
    install_requires=["numpy>=1.8", "scipy>=0.10",
                      "pymatgen>=3.1.5", "fireworks>=1.07",
                      "custodian>=0.8.2", "pymatgen-db>=0.4.8",
                      "ase>=3.9.0"],
    extras_require={"plot": ["matplotlib>=1.4.2"],
                    "babel": ["openbabel", "pybel"],
                    "remote": ["fabric"],
                    "doc": ["sphinx>=1.3.1", "sphinx-rtd-theme>=0.1.8"]
                   },    
    author = "Kiran Mathew, Joshua Gabriel, Arunima Singh, Richard G. Hennig",
    author_email = "km468@cornell.edu",
    description = ("High throughput analysis of interfaces using VASP and Materials Project tools"),
    license = "GPL",
    url = "https://github.com/henniggroup/MPInterfaces",
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
