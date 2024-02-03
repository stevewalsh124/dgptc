#!/bin/bash
#SBATCH -N1 --ntasks-per-node=4
#SBATCH -t 24:00:00
#SBATCH -p normal_q
#SBATCH -A Precipit

module purge
module load gcc/8.2.0
module load apps site/tinkercliffs/easybuild/setup
module load R/4.1.0-foss-2021a

#export R_LIBS="$HOME/R/gcc/3.6/" # dir for gcc version of package

export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

i=${1:-1}
#j=${2:-0}

R CMD BATCH --no-restore "--args mte=$i ncores=$OMP_NUM_THREADS" 1D_real_study_full_emuspace.R Rout/1D_realemufull${i}.Rout
