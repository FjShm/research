#!/bin/bash

./Crystallinity-A.o param.yaml
echo crystallinity-A.o was done.
echo -e "0.1\n5" | python plot.py sample_data/cry.txt
