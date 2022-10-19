#!/bin/bash

./Scattering.o param.yaml
python plot.py sample_data/out.scatter.txt
