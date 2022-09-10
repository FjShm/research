#!/bin/bash

./Crystallinity-A.o param.yaml
echo crystallinity-A.o was done.
python plot.py sample_data/cry.txt
