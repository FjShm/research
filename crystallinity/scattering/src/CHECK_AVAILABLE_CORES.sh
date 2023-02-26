#!/bin/bash -eu

function ql (){
    qstat -f | grep BIP | grep -v au | grep -v E | awk "{cnt=split(\$3,cores,\"/\"); s=cores[3]-cores[2]; if(s!=0){printf \"%20s\t%2d\n\",\$1,s}}" | sort -r -n -k 2,2
}

if [ $# -ne 1 ]; then
    echo "Usage: ./CHECK_AVAILABLE_CORES.sh [node-name]"
    exit 1
fi

echo -n "TOTAL CORES: "
dat=`ql | grep $1`
awk '{a+=$2} END{print a;}' <<< "$dat"
echo ---------------------
echo "$dat"

exit 0
