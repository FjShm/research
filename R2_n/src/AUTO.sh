#!/bin/bash

./R2n.o param.yaml
echo "R2n.o finished."
python plot.py sample_data/r2_n.txt
