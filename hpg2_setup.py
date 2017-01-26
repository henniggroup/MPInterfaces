import glob
import os

from setuptools import setup, find_packages

MPINT_DIR = os.path.dirname(os.path.abspath(__file__))

setup(
    name="mpinterfaces",
    version="1.2.0",

    install_requires=["FireWorks==1.3.3",
                      "custodian==1.0.1", "pymatgen-db==0.5.1",
                      "ase==3.11.0", "six"],
    extras_require={"babel": ["openbabel", "pybel"],
                    "remote": ["fabric"],
                    "doc": ["sphinx>=1.3.1", "sphinx-rtd-theme>=0.1.8"]
                    },
    author="Kiran Mathew, Joshua Gabriel, Arunima Singh, Richard G. Hennig",
    author_email="km468@cornell.edu",
    description=(
        "High throughput analysis of interfaces using VASP and Materials Project tools"),
    license="MIT",
    url="https://github.com/henniggroup/MPInterfaces",
    packages=find_packages(),
    long_description=open(
        os.path.join(os.path.dirname(__file__), 'README.rst')).read(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    scripts=glob.glob(os.path.join(MPINT_DIR, "scripts", "*"))
)
