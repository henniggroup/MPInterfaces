**Setup instructions on UF HPG2**

* Download and install the latest virtualenv at  https://pypi.python.org/pypi/virtualenv#downloads

* Check that you have python/2.7.8 loaded in your modules

* In your directory of choice under your home, install the virtualenv with
  python <path_to_virtualenv_extracted_folder>/virtualenv.py <path_to_where_you_want_virtualenv_with_name>.

  Lets call it mpint_venv

* Activate it with
  source mpint_venv/bin/activate

* pip install --ignore  numpy, monty

Finally,

* python setup.py develop


If this does not work please contact me at joshgabriel92@gmail.com

---- Try installing mpinterfaces on python3 venv at a version that works and syncs with matplotlib
      - use setup_hpg2.py script fot this
      - 1 . module load python3/3.5            ## switch to python3 module on hpg2
      - 2.  git clone https://github.com/materialsproject/pymatgen.git   ## select their most stable release
       which is v4.6.0
      - 3. virtualenv -p python3 <venv_name>
      - 4. cd mpinterfaces
      - 5. python setup_hpg2.py develop

Anaconda-Docker install only 1 issue:

-- how to ship the potcar files and/or point to them ?  
