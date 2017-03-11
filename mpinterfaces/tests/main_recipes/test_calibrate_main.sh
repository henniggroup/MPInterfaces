#!/bin/bash

# If there exists a configuration file already, store it.
cp ../../config_mine.yaml config_mine_saved.yaml

cp ../../../config.yaml ../../config_mine.yaml

python ../../calibrate.py > Calibrate_py_output.txt

File_differs_from_test(){
    if [[ $(diff $1 \
          ../../test_files/$1 2>&1|\
           tr -d ' \n\r\t ' | wc -c) -ne 0 ]]; then
        echo "Error: $1 does not match test file or test file does not exist."
        exit 1
    fi
    rm $1
}

File_differs_from_test "Calibrate_py_output.txt"
#File_differs_from_test "calibrate.json"

mv config_mine_saved.yaml  ../../config_mine.yaml 
rm -rf Slab* custodian.json
