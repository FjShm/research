#!/bin/bash

if [ ! -d results/data ]; then
    echo "simulation is not submitted."
    exit 1
fi

cd results/data
if [ -d joined ]; then
    echo "'results/data/joined' exist."
    echo "excute 'rm -r results/data/joined' first."
    exit 1
fi
mkdir joined
cd 1
logfilenames=(`ls | grep "LOG.*.out"`)
cd ..

for f in ${logfilenames[@]}
do
    echo joining $f...
    find . -name "$f"
    paste `find . -name "$f"` -d "" > .tmp
    tac .tmp > pmt.
    cat .tmp pmt. > joined/$f
    rm .tmp pmt.
done

exit 0
