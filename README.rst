About:
========

Python package that enables high throughput analysis of interfaces(=surface+ligand+solvent) using VASP_, VASPsol_ and materialsproject_ tools

.. _materialsproject: https://github.com/materialsproject

.. _VASPsol: http://vaspsol.mse.ufl.edu/

.. _VASP: http://www.vasp.at/

Packages:
==========

- mpinterfaces

Scripts:
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

- get the latest version from bitbucket
  
  if starting from scratch:
	
  * git clone https://username@bitbucket.org/matk86/vasp_automation.git

  else:

  * cd vasp_automation

  * git pull
	
- cd vasp_automation
	
- python setup.py install(or develop)

Status:
=======================

	Alpha testing phase.

Authors
=======================
   
	Kiran Mathew
	
	Joshua Gabriel
	
