#!/bin/bash

./Crystallinity-B.o param.yaml
echo crystallinity-B.o was done.
python plot.py sample_data/cry.txt
