#!/bin/bash

./MSD.o param.yaml
if [ $? -eq 0 ]; then
    echo "MSD.o finished."
else
    echo "MSD.o failed..."
fi
