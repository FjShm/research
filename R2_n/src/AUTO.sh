#!/bin/bash

./R2n.o param.yaml
if [ $? -eq 0 ]; then
    echo "R2n.o finished."
    python plot.py sample_data/r2_n.txt
else
    echo "R2n.o failed..."
fi
