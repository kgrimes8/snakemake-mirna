#!/bin/bash

if { conda env list | grep -i "testenv"; } &> /dev/null
then
    echo "conda testenv exists already"
else
    conda env create -f testenv.yaml
fi

conda run -n testenv pytest test-scripts.py