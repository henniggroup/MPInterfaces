#!/bin/bash

if [[ $(diff ../../POSCAR_diacetate_boxed.vasp \
      ../../test_files/POSCAR_diacetate_boxed.vasp |\
       tr -d ' \n\r\t ' | wc -c) -ne 0 ]]; then
    echo "Error: POSCAR_diacetate_boxed.vasp does not match test file"
    exit 1
fi
    echo "Files match: success"
