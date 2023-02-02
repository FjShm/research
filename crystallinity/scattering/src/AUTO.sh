#!/bin/bash

# input
DUMP_PATH=../../../sample_data/erate_0.1_dt_5.u.lammpstrj
ROTATIONTXT_PATH=../../../sample_data/rotation.txt
ODIR="."
ASPECT=xz

DUMP_FRAMES=6
N=49
M=512

RATIO=(0 0.5 1.0)

CORES=4
K=2
RESOLUTION=101

NODE=ika1.q

# functions
function round (){
    if [ $# -ne 2 ]; then
        echo "round [float] [digit]"
        exit
    fi

    root=`bc -l <<< 10^$2`
    if [ `bc -l <<< "$1 < 0"` -eq 1 ]; then
        number_tmp=`bc -l <<< $1*$root-0.5`
    else
        number_tmp=`bc -l <<< $1*$root+0.5`
    fi
    number_tmp=`awk "BEGIN {printf(\"%d\", $number_tmp)}"`
    #number_tmp=`bc -l <<< $number_tmp/$root`
    number_tmp=`awk "BEGIN {printf(\"%f\", $number_tmp/$root)}"`
    echo $number_tmp
}

function linspace (){
    # similar to numpy.linspace
    if [ $# -ne 4 ]; then
        echo "linspace [min] [max] [num] [dtype]"
        exit
    fi

    _dt_=`bc -l <<< "($2-($1))/($3-1)"`
    for _i_ in `seq 0 $(($3-1))`
    do
        val=`bc -l <<< $1+$_dt_*$_i_`
        if [ $4 == "float" ]; then
            val=`round $val 7`
            echo "$val"
        elif [ $4 == "int" ]; then
            val=`round $val 0`
            echo "${val%.*}"
        fi
    done
}

function reshape_for_ky (){
    # 降順で, 0以上の数を抽出
    _arr_=($1)
    _num_=${#_arr_[@]}
    for _i_ in `seq $(($_num_-1)) -1 0`
    do
        if [ `bc -l <<< "${_arr_[$_i_]} >= 0"` -eq 1 ]; then
            echo ${_arr_[$_i_]}
        fi
    done
}

function array_to_yamllist (){
    _arr_=($1)
    for str in ${_arr_[@]}
    do
        echo -n "\n  - $str"
    done
}


# prepare
if [ $CORES -gt $RESOLUTION ]; then
    echo "CORES has to be smaller than RESOLUTION"
    echo "fix CORES from $CORES to $RESOLUTION"
    CORES=$RESOLUTION
fi
if [ -d results/data ]; then
    echo "directory 'results/data' exists."
    echo "excute \`rm -r results/data\`"
    exit 1
fi
cwd=`pwd`
mkdir -p results/data && cd results/data
kx_all=(`linspace -$K $K $RESOLUTION float`)
ky=(`reshape_for_ky "${kx_all[*]}"`)
## core-n has jobs of idx[n] <= id < idx[n+1]
idx=(`linspace 0 $RESOLUTION $(($CORES+1)) int`)
RATIO=`array_to_yamllist "${RATIO[*]}"`
ky=`array_to_yamllist "${ky[*]}"`

# make dir and copy files
for c in `seq 1 $CORES`
do
    echo -e "\rnow $c/$CORES..."
    mkdir $c
    cp $cwd/Scattering.o $cwd/run.sh $c
    cp $cwd/format_param.yaml $c/param.yaml
    cd $c
    kx_idx=(`seq ${idx[$(($c-1))]} $((${idx[$c]}-1))`)
    kx_idx=`array_to_yamllist "${kx_idx[*]}"`
    kx=(`echo "${kx_all[@]: ${idx[$(($c-1))]}: ${idx[$c]}}"`)
    kx=`array_to_yamllist "${kx[*]}"`
    # overwrite
    sed -i -e "s;__DUMP_PATH__;$DUMP_PATH;g" param.yaml
    sed -i -e "s;__ROTATIONTXT_PATH__;$ROTATIONTXT_PATH;g" param.yaml
    sed -i -e "s;__ODIR__;$ODIR;g" param.yaml
    sed -i -e "s;__ASPECT__;$ASPECT;g" param.yaml
    sed -i -e "s;__DUMP_FRAMES__;$DUMP_FRAMES;g" param.yaml
    sed -i -e "s;__N__;$N;g" param.yaml
    sed -i -e "s;__M__;$M;g" param.yaml
    sed -i -e "s;__KX__;$kx;g" param.yaml
    sed -i -e "s;__KY__;$ky;g" param.yaml
    sed -i -e "s;__IDX__;$kx_idx;g" param.yaml
    sed -i -e "s;__RATIO__;$RATIO;g" param.yaml
    sed -i -e "s;__NODE__;$NODE;g" run.sh
    cd ..
done


exit 0
