import argparse
from numpy import arange
import os
import pathlib
import time
import sys
import subprocess
import shutil

def to_list(arg):
    return [float(i) for i in arg.split(",")]


def to_string_list(arg):
    return [str(i) for i in arg.split(",")]

# doc: save run information to run log file, append if necessary
def saveLogs(cmdinput) -> str:
    options = f"start time: {time.ctime()}\n"
    options += "  + command: python " + f"{' '.join(sys.argv)}\n"
    options += "\n".join("  + {}: {}".format(k, v) for k, v in cmdinput.items())
    print(options)

    savepath = cmdinput["data_path"]
    pathlib.Path(f"{savepath}/logs").mkdir(parents=True, exist_ok=True)
    simInfoFile = f"{savepath}/logs/search_table_information.log"
    with open(simInfoFile, mode="a") as historyFile:
        historyFile.write(options)
    return simInfoFile

parser = argparse.ArgumentParser()
parser.add_argument("--numberTrials", type=int, default=1, help="number of trials to run")
parser.add_argument("--numberSamples", type=int, default=50, help="number of samples to generate")
parser.add_argument("--dims", required=True, type=to_list, help="dimensions of the lattice")
parser.add_argument("--data_path", required=True, help="path to save data")

parser.add_argument("--intensity", type=float, default=-1, help="hotspot intensity")
parser.add_argument("--radius", type=int, default=-1, help="hotspot radius")
parser.add_argument("--density", type=float, default=-1, help="hotspot density")

parser.add_argument("--gap", type=int, default=0, help="empty space around the lattice edges")
parser.add_argument("--ref_line", type=int, default=0, help="to remove lattice effects at the front, one may desire to terminate data colection at a certain line")

parser.add_argument("--rng_seed", type=int, default=1, help="random seed")
parser.add_argument("--detailed_analytics", action="store_true", help="save lineages, sectors, branch points, and final front")
parser.add_argument("--rewrite", action="store_true", help="overwrite existing data during generation")

parser.add_argument("--nEnvs", type=int, default=1, help="number of environments to run")
parser.add_argument("--parameters", default="density,intensity", type=to_string_list, help="the two parameters to vary. Choose from ['density', 'intensity', 'radius']. The non-chosen parameter will be kept constant, and must be specified.")
# task: append multiple intervals when calling flags repeatedly
parser.add_argument("--intervals_1", required=True, type=to_list, help="intervals for the first parameter")
parser.add_argument("--intervals_2", required=True, type=to_list, help="intervals for the second parameter")

args = parser.parse_args()

extra_flags = []
if args.detailed_analytics:
    extra_flags.append("detailed_analytics")
if args.rewrite:
    extra_flags.append("rewrite")

p_1, p_2 = args.parameters
start, step, stop = args.intervals_1

intensity = args.intensity
radius = args.radius
density = args.density

parameter_options = ["intensity", "radius", "density"]
third_parameter = [x for x in parameter_options if x not in [p_1, p_2]]

print(f"Varying {p_1} and {p_2} over {args.intervals_1} and {args.intervals_2}, respectively")
print(f"Keeping {third_parameter[0]} constant")

if vars(args)[third_parameter[0]] < 0:
    print(f"Please specify the value of {third_parameter[0]}")
    exit()

#simFileInfo = saveLogs(vars(args))

for parameterVal in arange(start, stop + step, step=step):
    if p_1 == "intensity":
        intensity = parameterVal
    elif p_1 == "radius":
        radius = parameterVal
    elif p_1 == "density":
        density = parameterVal
    else:
        print("unable to parse your commands. Examine script/cmd")
        exit()

    cmd = f"""sbatch src/processing/slurm_main.sh \
            --data_path {args.data_path} \
            --numberTrials {args.numberTrials} \
            --numberSamples {args.numberSamples} \
            --dims {",".join(str(x) for x in args.dims)} \
            --intensity {intensity} \
            --density {density} \
            --radius {radius} \
            --parameter {p_2} \
            --nEnvs {args.nEnvs} \
            --gap {args.gap} \
            --ref_line {args.ref_line} \
            --intervals {",".join(str(x) for x in args.intervals_2)} \
            --flags {",".join(extra_flags)}
            """

    if shutil.which("slurm") is None:
        print("Slurm scheduler is not installed. Executing using bash.")
        cmd = cmd.replace("sbatch", "bash")

    subprocess.run(cmd, shell=True)
