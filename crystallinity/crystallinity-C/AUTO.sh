#!/bin/bash

./Crystallinity-C.o param.yaml
echo crystallinity-C.o was done.
python plot.py sample_data/cry.txt
