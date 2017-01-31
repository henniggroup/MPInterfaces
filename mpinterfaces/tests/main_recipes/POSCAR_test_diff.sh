#!/bin/bash

File_differs_from_test(){
    echo "Testing $1"
    if [[ $(diff ../../$1 \
          ../../test_files/$1 |\
           tr -d ' \n\r\t ' | wc -c) -ne 0 ]]; then
        echo "Error: $1 does not match test file"
        exit 1
    fi
        echo "Files match: success"
}

File_differs_from_test "POSCAR_diacetate_boxed.vasp"
File_differs_from_test "POSCAR_interface.vasp"
File_differs_from_test "POSCAR_slab.vasp"
