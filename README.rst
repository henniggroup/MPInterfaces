About:
========

Python package that enables high throughput analysis of interfaces(ligands/solvent/hetero-structures) using VASP_, VASPsol_ and materialsproject_ tools

.. _materialsproject: https://github.com/materialsproject

.. _VASPsol: http://vaspsol.mse.ufl.edu/

.. _VASP: http://www.vasp.at/

Packages:
==========

- mpinterfaces

Examples:
==========

- sample scripts that demonstrate the usage of mpinterfaces as well as materialsproject packages
- start from here

Requirements:
==============

- python>=2.7
- numpy
- scipy
- pymatgen
- custodian
- fireworks
- pymatgen-db

..

	It is highly recommended that you install this package and all the required packages in a virtual environment( http://virtualenv.readthedocs.org/en/latest/virtualenv.html)

Install Instructions:
=======================

- get the latest version from github
  
- cd MPInterfaces
	
- python setup.py install(or develop)

- to get the documentation (sphinx package must be available)

  * html format
    
    - cd docs; make html

  * pdf

    - cd docs; make latexpdf

Status:
=======================

	Alpha testing phase.

Authors
=======================
   
	Kiran Mathew
	
	Joshua Gabriel

	Richard G. Hennig
