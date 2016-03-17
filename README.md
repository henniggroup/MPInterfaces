# twod_materials
Python modules for high-throughput 2D Materials Characterization

# Setup
Edit config.yaml with your system's settings, following the template:

```
username: your_hipergator_username
mp_api: your_materials_project_api_key
normal_binary: path_to_normal_vasp_executable
twod_binary: path_to_twod_vasp_executable
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

Please ensure that a compiled version of vasp exists in your ~/bin, and
that a specially compiled vasp with no relaxation in the z-direction also
exists in your ~/bin. They should be called `vasp` and `vasp_noz`,
respectively.

See the examples folder for some sample workflows.
