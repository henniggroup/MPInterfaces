## Overview of Tests ##

The following modules MUST NOT break with any additional feature updates
or API improvements

* __calibrate.py__ writes out a list of directories, when python calibrate.py is
  executed. This is essential to the creation and submission of jobs to the
  queue.

  A passed test is one where all the functions called in main work correctly
  together and is a more rigorous pass than individual functions passing
  the test

  Running this test requires having access to copyrighted POTCAR files and
  using your config_mine.yaml to state where they are. Therefore, this test is
  skipped on all automated testing.

* __interface.py__ writes out poscar files specified in the test_files
  directory. There should be a complete match of the files.

  This is again a test that passes if the complete call of functions as
  defined by main passes

* __transformations.py__ is primarily a group of functions used to create the
  heterostructure interface. A successful pass is a successful run of the
  __heterostructure.py__ examples script

* a complete __package-production-pass__ is when __examples/vasp_jobs.py__ successfully runs
  and submits jobs to the queue and checkpoint files get updated correctly for
  a test workflow

  the successful execution of __examples/vasp_jobs.py__ means a successful execution
  of the helping functions in the modules calibrate.py, instrument.py and utils.py
