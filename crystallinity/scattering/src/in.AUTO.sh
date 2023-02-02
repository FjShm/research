#!/bin/bash -eu

DUMP_PATH=../../../sample_data/erate_0.1_dt_5.u.lammpstrj
ROTATIONTXT_PATH=../../../sample_data/rotation.txt
ODIR="."
ASPECT=xz

DUMP_FRAMES=6
N=49
M=512

RATIO=(0 0.5 1.0)

K=2
RESOLUTION=101

declare -A NODES=(
    ["ika1.q"]=4
    ["tuna0.q"]=4
)
