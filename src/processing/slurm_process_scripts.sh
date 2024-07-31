#!/bin/bash
#SBATCH --job-name=metrics
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=10gb
#SBATCH --time=336:00:00
#SBATCH --output=workspace/logs/slurm/metrics_%j.log
pwd; hostname; date

pararser() {
    # default values
    ncores=${ncores:-"3"}
    table=${table:-"workspace/sim_collections/parameter_space_table.csv"}
    script=${script:-"src/processing/process_lineageMSD.jl"}
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

# get opts options
pararser $@

IFS="," read -ra ADDR <<< "$flags"
for flg in "${ADDR[@]}"; do 
    extra_flags+=" --$flg"
done

echo $extra_flags

echo "Running your script now"

echo "--script $script --parameter_space_table $table $extra_flags" 

julia --project=@. -O3 -t $ncores $script --parameter_space_table $table $extra_flags

echo "Completed running your script."
date
