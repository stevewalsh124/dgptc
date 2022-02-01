#!/bin/bash

# Make sure to change the -A line to your allocation

#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -t 48:00:00
#SBATCH -p normal_q
#SBATCH -A Precipit

# Add modules
module purge
module load singularity

# Change to the directory from which the job was submitted
#cd $PBS_O_WORKDIR

# Run R


i=${1:-1}
#j=${2:-0}

#nohup R CMD BATCH "--args seed=$i reps=$reps" spam_mc.R spam_mc_$i.Rout &

echo "$( date ): Starting TC $i now"
singularity exec --bind=/groups:/groups /groups/arcsingularity/ood-rstudio-geospatial_4.0.3.sif R CMD BATCH --no-restore "--args ste=$i" /home/walsh124/dgptc/2D_examples/2D_TCexample.R Rout/prediction${i}.Rout 
echo "$( date ): Finished prediction $i" 

exit
