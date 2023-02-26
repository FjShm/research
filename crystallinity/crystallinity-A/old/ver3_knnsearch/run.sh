#!/bin/sh
#$ -S /bin/sh
#$ -N aris
#$ -cwd
#$ -j y
#$ -q ika1.q
#$ -pe openmp 2


./AUTO.sh
