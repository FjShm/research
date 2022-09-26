#!/bin/sh
#$ -S /bin/sh
#$ -N aris
#$ -cwd
#$ -j y
#$ -q ara.q
#$ -pe openmp 1


./Show_crystallinity.o param.yaml
