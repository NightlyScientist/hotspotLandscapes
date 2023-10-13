#!/bin/bash
#SBATCH --job-name=parameter_search
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=discardmyinformation@gmail.com 
#SBATCH --ntasks=10
#SBATCH --mem=20gb
#SBATCH --time=48:00:00
#SBATCH --output=/home/jgonzaleznunez/Projects/Research/Main/lineage-landscapes/storage/scratch/sbatch_logs/submission_%j.log   # Standard output and error log

pwd; hostname; date

pararser() {
    # default values
    radius=${radius:-"10"}
    density=${density:-"0.1"}
    intensity=${intensity:-"0"}
    parameter=${parameter:-"density"}
    intervals=${intervals:-"0,0.1,1"}
    dims=${dims:-"2000,1000"}
    numberTrials=${numberTrials:-"1200"}
    nCollections=${nCollections:-"2"}

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
logFile="/home/jgonzaleznunez/Projects/Research/Main/lineage-landscapes/storage/scratch/logs/submissions.log"

# get cli options
pararser $@

echo "Running your script now"
echo

echo "radius:$radius | density: $density | intensity: $intensity"
echo "parameter:$parameter | intervals: $intervals"
echo "numberTrials: $numberTrials | nCollections: $nCollections"

echo

# add sumbission parameters to logs
echo "start-------------" >> $logFile
echo "radius: $radius | density: $density | intensity: $intensity" >> $logFile
echo "parameter: $parameter | intervals: $intervals" >> $logFile
echo "numberTrials: $numberTrials | nCollections: $nCollections" >> $logFile
echo "dimensions: $dims" >> $logFile
echo "end--------------" >> $logFile

python main.py --numberTrials $numberTrials --numberSamples 400 --dims $dims --density $density --intensity $intensity --radius $radius --savepath storage/data/production/ --parameter $parameter --intervals=$intervals --nCollections $nCollections --num_threads 10 --background

date
