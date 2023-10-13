#!/bin/bash
#SBATCH --job-name=msd
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=thelearningprofile@gmail.com
#SBATCH --ntasks=2
#SBATCH --mem=20gb
#SBATCH --time=72:00:00
#SBATCH --array=0-6%4
#SBATCH --output=/home/jgonzaleznunez/Projects/Research/Main/lineage-landscapes/storage/scratch/sbatch_logs/submission_%j.log
pwd; hostname; date

echo "Running your script now"

# bash array of all directories
DIRS=($1/*) 

echo ${DIRS[$SLURM_ARRAY_TASK_ID]} 

julia --project=@. src/analysis/lineageMSD/lineageMSD.jl --path ${DIRS[$SLURM_ARRAY_TASK_ID]}  

date