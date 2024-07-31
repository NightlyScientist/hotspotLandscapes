#!/bin/bash
#SBATCH --job-name=mut
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00
#SBATCH --output=workspace/logs/slurm/submission_%j.txt

pwd; hostname; date

pararser() {
    # default values
    density=${density:-"0.07"}
    radius=${radius:-"10"}
    intensity=${intensity:-"0.0"}
    
    parameter=${parameter:-"intensity"}
    intervals=${intervals:-"0,0.1,1"}
    
    dims=${dims:-"500,1000"}
    numberTrials=${numberTrials:-"50"}
    numberSamples=${numberSamples:-"50"}

    ref_line=${ref_line:-"0"}
    gap=${gap:-"0"}
    nEnvs=${nEnvs:-"1"}

    data_path=${data_path:-"workspace/sims"}
    flags=${flags:-""}

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
logFile="workspace/logs/jobs/submission_record.log"
touch $logFile

# get cli options
pararser $@

echo "Running your script now"
echo

echo "parameter:$parameter | intervals: $intervals"
echo "numberTrials: $numberTrials"

echo

# add sumbission parameters to logs
echo "start-------------" >> $logFile
echo "parameter: $parameter | intervals: $intervals" >> $logFile
echo "numberTrials: $numberTrials" >> $logFile
echo "dimensions: $dims" >> $logFile
echo "end--------------" >> $logFile
echo "       " >> $logFile

# parse extra flags that evaluate to 'store_true'
IFS=',' read -ra ADDR <<< "$flags"
for flg in "${ADDR[@]}"; do
  extra_flags+=" --$flg" 
done

echo $extra_flags

python src/processing/main.py --numberTrials $numberTrials --numberSamples $numberSamples --dims $dims --intensity $intensity --radius $radius --density $density --data_path $data_path --parameter $parameter --intervals=$intervals --num_threads 4 --nEnvs $nEnvs --ref_line $ref_line --gap $gap --background $extra_flags


date
