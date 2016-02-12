# twod_materials
Python modules for high-throughput 2D Materials Characterization

# Setup
Edit config.yaml with your system's settings, following the template:

```
mp_api: your_materials_project_api_key
potentials: path_to_your_vasp_potentials
```

And then move it to your home directory ```mv config.yaml```

# Usage

This package is designed to characterize 2D materials with as little 
*a priori* knowledge as possible. The only input required are the 
structures of the materials, which should be stored as POSCAR files in 
uniquely named directories, e.g.:

+ My_2D_Search
    + WSe2
        + POSCAR
    + MoS2
        + POSCAR
    + Ti2CO2
        + POSCAR

In almost all cases, before doing anything else, the user should use the 
relax() function to optimize the structures of the 2D materials using 
the framework of input parameters included here.

See the examples folder for some sample workflows.