#!/bin/bash -eu

# used in param.yaml
## dump path from me
DUMP_PATH=sample_data/erate_0.1_dt_5.u.lammpstrj

## rotation.txt path from me
ROTATIONTXT_PATH=sample_data/rotation.txt

ODIR="."

ASPECT=xz

DUMP_FRAMES=6
FRAME_NUMBERS=(1 3 6)
N=49
M=512


K=2
RESOLUTION=501


# used in qsub.sh
declare -A NODES=(
    ["ika1.q"]=60
    ["tuna0.q"]=120
)

# used in AUTO.sh
HEADER=9
