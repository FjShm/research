#!/bin/bash

./Scattering.o param.yaml
echo "Scattering.o was finished!!"
/home/fujii/anaconda3/bin/python3 plot.py sample_data/out.scatter.txt
