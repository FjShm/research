#!/bin/bash

./Scattering.o param.yaml

if [ $? -eq 0 ]; then
    echo "Scattering.o was finished!!"
    /home/fujii/anaconda3/bin/python3 plot.py sample_data/out.scatter.txt
else
    echo "Scattering.o failed..."
fi
