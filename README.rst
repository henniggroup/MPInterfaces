Python package that uses the materialsproject tools to automate vasp calculations for molecules, surfaces and interfaces(=surface+ligand+solvent)

Packages:
==========

- mpinterfaces

scripts:
==========

- some utility scripts that make use of the modules in the package

Requirements:
==============

- python2.7
- numpy
- scipy
- pymatgen
- custodian
- fireworks
- pymatgen-db

..

	It is highly recommended that you install this package and all the required packages in a virtual environment

Install Instructions:
=======================

- fetch the package from bitbucket
  
  if starting from scratch:
	
  * git clone https://username@bitbucket.org/matk86/vasp_automation.git

  else:

  * git pull
	
- cd vasp_automation
	
- python setup.py install (no sudo privilege required if installing in a virtual environment)

Status:
=======================

	Alpha testing phase. Requires more thorough testing

:Authors:
	Kiran Mathew
	Joshua Gabriel
	
