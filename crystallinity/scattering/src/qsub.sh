#!/bin/sh
#$ -S /bin/sh
#$ -N aris
#$ -cwd
#$ -j y
#$ -q __NODE__
#$ -pe openmp 1


./Scattering.o param.yaml
