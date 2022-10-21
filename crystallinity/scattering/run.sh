#!/bin/sh
#$ -S /bin/sh
#$ -N aris
#$ -cwd
#$ -j y
#$ -q sujiko.q
#$ -pe openmp 52
##$ -o out.txt

#$ -v KMP_AFFINITY=verbose,compact,1
#$ -v MKL_NUM_THREADS=104
#$ -v OMP_NUM_THREADS=104
#$ -v MKL_DOMAIN_NUM_THREADS="MKL_ALL=1, MKL_FFT=1"
#$ -v MKL_DYNAMIC=FALSE
#$ -v OMP_DYNAMIC=FALSE
#$ -v OMP_SCHEDULE="static"

echo "=== START ==="
date

#export I_MPI_PIN=0
export I_MPI_DEBUG=1
export FI_PROVIDER=tcp
export I_MPI_PIN=1
export I_MPI_DEVICE=ssm
export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_MPD_RSH=/usr/bin/rsh
export I_MPI_MPD_TMPDIR=/tmp
source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1

./AUTO.sh

echo "=== END ==="
date
