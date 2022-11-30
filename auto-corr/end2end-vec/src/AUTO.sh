#!/bin/bash

./Ete_auto_corr.o param.yaml
if [ $? -eq 0 ]; then
    echo "Ete_auto_corr.o finished."
else
    echo "Ete_auto_corr.o failed..."
fi
