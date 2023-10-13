#!/bin/bash
#SBATCH --job-name=shortestPaths
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=thelearningprofile@gmail.com
#SBATCH --ntasks=2
#SBATCH --mem=20gb
#SBATCH --time=18:00:00
#SBATCH --output=/home/jgonzaleznunez/Projects/Research/Main/lineage-landscapes/storage/scratch/sbatch_logs/submission_%j.log

pararser() {
    # default values
    input=${input:-""}
    output=${output:-"storage/data/predictions/"}

    # Assign the values given by the user
    while [ $# -gt 0 ]; do
        if [[ $1 == *"--"* ]]; then
            param="${1/--/}"
            declare -g $param="$2"
        fi
        shift
    done
}
# where to save log information
logFile="/home/jgonzaleznunez/Projects/Research/Main/lineage-landscapes/storage/scratch/logs/predictions.log"

# get cli options
pararser $@

pwd; hostname; date

echo "Running your script now"
echo

echo "input: $input"
echo "output: $output"

echo

echo "start-------------" >> $logFile
echo "input: $input" >> $logFile
echo "output: $output" >> $logFile
echo "end--------------" >> $logFile

julia --threads 2 --project=@. src/predictions/shortestpaths.jl  --output $output --input $input

date
