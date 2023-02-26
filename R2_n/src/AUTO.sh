#!/bin/bash

./R2n.o param.yaml
if [ $? -eq 0 ]; then
    echo "R2n.o finished."
    ~/anaconda3/bin/python plot.py sample_data/R2_n.dat
else
    echo "R2n.o failed..."
fi
