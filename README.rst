MPinterfaces is a python package that enables high throughput Density
Functional Theory(DFT) analysis of arbitrary material interfaces(ligand capped
nanoparticles, surfaces in the presence of solvents and hetero-structure
interfaces) using VASP_, VASPsol_ and materialsproject_ database as well
as their open source tools.

.. _materialsproject: https://github.com/materialsproject

.. _VASPsol: https://github.com/henniggroup/VASPsol

.. _VASP: http://www.vasp.at/


Requirements
==============

1. Make sure that you are using python>=2.7 (do a "python --version")

2. It is highly recommended that you use gcc compiler. So type::

   export CC=gcc

3. Unless you have admin privilege on the machine you are installing, it is
   better to install this package and all its dependencies in a virtual environment.

   - get the latest version from https://pypi.python.org/pypi/virtualenv#downloads
   
   - tar xvfz virtualenv-X.X.X.tar.gz
   
   - cd virtualenv-X.X.X
   
   - setup the virtual environment in ~/myvenv (or set to some other path and folder name)
     
     * python virtualenv.py ~/myvenv
   
   For detailed instructions and documentaion see
   http://virtualenv.readthedocs.org/en/latest/installation.html

4. Install numpy::

   pip install numpy

   
Installation
==============

Using pip
----------

pip install mpinterfaces


Bleeding edge
-------------

- Clone the latest version from github
  
- cd MPInterfaces
	
- python setup.py install(or develop)

  
Documentation
==============

A very minimal documentation is avaiable at the moment and work in underway
to improve it. We use the sphinx package to generate the documentation.
First install the package 'sphinx' and the theme package 'sphinx-rtd-theme'
using pip and then to generate the documentation either in html or pdf format
do the following:

  * html format
    
    - cd docs; make html

    - open _build/html/index.html in a webbrowser to see the documentation

  * pdf

    - cd docs; make latexpdf

    - The pdf will be in the _build/latex folder 

      
Examples
==========

The examples folder contain some sample scripts that demonstrate the
usage of mpinterfaces as well as materialsproject packages


Authors
=========
   
Kiran Mathew
	
Joshua Gabriel

Richard G. Hennig

