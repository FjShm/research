#!/bin/sh
#$ -S /bin/sh
#$ -N aris
#$ -cwd
#$ -j y
#$ -q tuna0.q
#$ -pe openmp 1


./Scattering.o param.yaml
