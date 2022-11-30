#!/bin/bash

./Rouse_modep.o param.yaml
if [ $? -eq 0 ]; then
    echo "Rouse_modep.o finished."
else
    echo "Rouse_modep.o failed..."
fi
