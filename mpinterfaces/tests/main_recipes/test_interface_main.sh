#!/bin/bash

# If there exists a configuration file already, store it.
cp ../../config_mine.yaml config_mine_saved.yaml

cp ../../../config.yaml ../../config_mine.yaml

python ../../interface.py > Interface_py_output.txt

File_differs_from_test(){
    if [[ $(diff $1 \
          ../../test_files/$1 2>&1|\
           tr -d ' \n\r\t ' | wc -c) -ne 0 ]]; then
        echo "Error: $1 does not match test file or test file does not exist."
        exit 1
    fi
    rm $1
}

File_differs_from_test "Interface_py_output.txt"
File_differs_from_test "POSCAR_diacetate_boxed.vasp"
File_differs_from_test "POSCAR_interface.vasp"
File_differs_from_test "POSCAR_slab.vasp"
File_differs_from_test "lead_acetate.xyz"

mv config_mine_saved.yaml  ../../config_mine.yaml
