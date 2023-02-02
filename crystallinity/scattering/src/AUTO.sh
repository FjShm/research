#!/bin/bash -eu

SECONDS=0

# input
echo -n "[`date +\"%Y-%m-%d %H:%M:%S\"`] include $1..."
. $1
echo " done"

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

function monitor_job (){
    while :
    do
        qstat | grep $1 > /dev/null
        if [ $? -eq 1 ]; then
            break
        fi
        sleep 10s
    done
}


# validation
echo -n "[`date +\"%Y-%m-%d %H:%M:%S\"`] checking environments..."
## output directory
if [ -d results/data ]; then
    echo -e " failed\ndirectory 'results/data' exists."
    echo "excute \`rm -r results/data\`"
    exit 1
fi

## node names
for NODE in ${!NODES[@]}
do
    qstat -f | grep $NODE@ > /dev/null
    if [ $? -eq 1 ]; then
        echo Invalid node name: $NODE
        exit 1
    fi
done
echo " done"

# prepare
echo -n "[`date +\"%Y-%m-%d %H:%M:%S\"`] preparing some variables..."
cwd=`pwd`
mkdir -p results/data && cd results/data && mkdir joined
kx_all=(`linspace -$K $K $RESOLUTION float`)
ky=(`reshape_for_ky "${kx_all[*]}"`)

## cores
CORES=0
for CORE in ${NODES[@]}
do
    CORES=`bc <<< $CORES+$CORE`
    if [ $CORES -gt $RESOLUTION ]; then
        echo "CORES has to be smaller than RESOLUTION"
        echo "fix CORES from $CORES to $RESOLUTION"
        CORES=$RESOLUTION
    fi
done

## core-n has jobs of idx[n] <= id < idx[n+1]
idx=(`linspace 0 $RESOLUTION $(($CORES+1)) int`)
RATIO=`array_to_yamllist "${RATIO[*]}"`
ky=`array_to_yamllist "${ky[*]}"`
echo " done"

## save kx, ky
echo -n "[`date +\"%Y-%m-%d %H:%M:%S\"`] exporting 'results/data/kx.txt', and 'results/data/ky.txt'..."
echo -n "" > kx.txt
echo -n "" > ky.txt
for i in `seq 1 ${#kx_all[@]}`
do
    echo -n "${kx_all[$(($i-1))]} " >> kx.txt
    echo "${kx_all[-$i]}" >> ky.txt
done
ky_has_0=false
if [ `bc -l <<< "${ky[-1]} == 0.0"` -eq 1 ]; then
    ky_has_0=true
fi
echo " done"

# make dir and copy files
echo "[`date +\"%Y-%m-%d %H:%M:%S\"`] parallelizing..."
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

    ## 自nodeが$NODESに含まれていれば-1しておく
    hostnode="`qstat -f | grep \`hostname\``"
    hostnode=(${hostnode// / })
    if [ -n ${NODES[${hostnode[0]%@*}]} ]; then
        NODES[${hostnode[0]%@*}]=`bc <<< ${NODES[${hostnode[0]%@*}]}-1`
    fi

    ## select node
    for NODE in ${!NODES[@]}
    do
        if [ ${NODES[$NODE]} -eq 0 ]; then
            continue
        fi
        sed -i -e "s;__NODE__;$NODE;g" run.sh
        NODES[$NODE]=`bc <<< ${NODES[$NODE]}-1`
        break
    done

    ## 最後はlocalでrun, それ以外はqsub
    if [ $c -eq $CORES ]; then
        echo run at `hostname`
        ./Scattering.o param.yaml > /dev/null &
    else
        res=`qsub run.sh`
        res_splited=(${res// / })
        jobID=${res_splited[2]}
        monitor_job $jobID &
        echo $res
    fi
    cd ..
done
echo "[`date +\"%Y-%m-%d %H:%M:%S\"`] parallelizing... done"

# wait for all background processes
echo -n "[`date +\"%Y-%m-%d %H:%M:%S\"`] wait for all processes..."
wait
echo " done"

# join all data
## search LOG.*.out
echo "[`date +\"%Y-%m-%d %H:%M:%S\"`] start to join data..."
cd 1
logfilenames=(`ls | grep "LOG.*.out"`)
cd ..

## join
for f in ${logfilenames[@]}
do
    echo joining $f...
    find . -name "$f"
    paste `find . -name "$f"` -d "" > .tmp
    tac .tmp > pmt.
    if [ $ky_has_0 == "true" ]; then
        mv pmt. pmt.bef
        tail -n +2 pmt.bef > pmt. # 1行目だけskipして出力
        rm pmt.bef
    fi
    cat .tmp pmt. > joined/$f
    rm .tmp pmt.
done
echo -e "[`date +\"%Y-%m-%d %H:%M:%S\"`] All done!\n\n"

i=$SECONDS
set +e
((sec=i%60, min=(i%3600)/60, hrs=i/3600))
set -e
timestamp=$(printf "%d:%02d:%02d" "$hrs" "$min" "$sec")
echo "-------------------------------------------------"
echo "Wall Time:  $timestamp"

exit 0
