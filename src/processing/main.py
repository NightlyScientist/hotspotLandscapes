import argparse
from numpy import arange
import pathlib
import os
import time
import multiprocessing as mp
import sys
import copy
from random import seed, randint


class KeyboardInterruptError(Exception):
    pass


class Parameters:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def collectParameters(opts) ->"list[Parameters]":
    parameterList = []
    start, step, stop = opts["intervals"]
    
    if opts["rng_seed"] < 1:
        # initialize rng with opts seed
        seed(opts["rng_seed"])

        # create multiple collections for each landscape config
        rng = randint(1, 1e10)
    else:
        rng = opts["rng_seed"]

    for env_num in range(0, opts["nEnvs"]):
        # step through range of parameter parameter values 
        for parameterVal in arange(start, stop + step, step=step):
            cdict = copy.deepcopy(opts)

            # modify scanned parameter and data_path
            cdict[opts["parameter"]] = parameterVal
            parameters = Parameters(**cdict)
            parameters.data_path = f'{cdict["data_path"]}/env_{env_num}'
            parameters.rng_seed = rng + env_num
            parameterList.append((parameters))
    return parameterList


def to_list(arg):
    return [float(i) for i in arg.split(",")]


def runJulia(prms: Parameters):
    cmd_rewrite = "--rewrite" if prms.rewrite else ""
    cmd_detailed_anal = "--detailed_analytics" if prms.detailed_analytics else ""

    command = f"""julia --project=@. -O3 src/base/main.jl \
            --numberTrials {prms.numberTrials} \
            --numberSamples {int(prms.numberSamples)} \
            --width {int(prms.dims[0])} \
            --height {int(prms.dims[1])} \
            --intensity {round(prms.intensity, 4)} \
            --radius {int(prms.radius)} \
            --density {float(prms.density)} \
            --rng_seed {int(prms.rng_seed)} \
            --ref_line {int(prms.ref_line)} \
            --gap {int(prms.gap)} \
            --data_path {prms.data_path} \
            --parameter {prms.parameter} \
            {cmd_rewrite} \
            {cmd_detailed_anal} \
            1> /dev/null
        """
    os.system(command)


def main(nworkers, params, simInfoFile):
    pool = mp.Pool(processes=nworkers)
    t1 = time.time()
    try:
        print("starting the pool map")
        pool.map(runJulia, params)
        pool.close()
        print("pool map complete")
    except KeyboardInterrupt:
        print("got ^C while pool mapping, terminating the pool")
        pool.terminate()
        print("pool is terminated")
    except Exception as e:
        print("got exception: %r, terminating the pool" % (e,))
        pool.terminate()
        print("pool is terminated")
    finally:
        print("joining pool processes")
        pool.join()
        print("join complete")

    msg = f"\nfinished. Elapsed wall time: {round((time.time()-t1)/3600, 3)} hour(s), at {time.ctime()}"
    print(msg)

    with open(simInfoFile, mode="a") as historyFile:
        historyFile.write(msg)

parser = argparse.ArgumentParser()
parser.add_argument("--numberTrials", type=int, default=1)
parser.add_argument("--numberSamples", type=int, default=50)

parser.add_argument("--dims", required=True, type=to_list)

parser.add_argument("--rng_seed", type=int, default=1)
parser.add_argument("--rewrite", action="store_true")
parser.add_argument("--data_path", required=True)

parser.add_argument("--detailed_analytics", action="store_true")

parser.add_argument("--intensity", required=True, type=float)
parser.add_argument("--radius", type=int, default=10)
parser.add_argument("--density", type=float, default=0.09)

parser.add_argument("--parameter", default="density")
parser.add_argument("--intervals", required=True, type=to_list)
parser.add_argument("--nEnvs", type=int, default=1)
parser.add_argument("--ref_line", type=int, default=0)
parser.add_argument("--gap", type=int, default=0)

parser.add_argument("--num_threads", type=int, default=10)
parser.add_argument("--background", action="store_true")

args = parser.parse_args()
opts = vars(args)

opts[args.parameter] = "xxx"
pc = Parameters(**opts)

names = [
    f"XY:{int(pc.dims[0])},{int(pc.dims[1])}",
    f'D:{round(pc.density, 3) if type(pc.density) != str else pc.density}',
    f'I:{round(pc.intensity, 2) if type(pc.intensity) != str else pc.intensity}',
    f'R:{round(pc.radius, 3) if type(pc.radius) != str else pc.radius}',
    f"nEnvs:{pc.nEnvs}",
    f'gap:{pc.gap}',
    f'intervals:{",".join([str(round(s, 3)) for s in pc.intervals])}',
    f"trials:{int(pc.numberTrials)}"]

if args.detailed_analytics:
    names.append("da")

opts[args.parameter] = "xxx"
pc = Parameters(**opts)

data_path = f"{args.data_path}/" + "_".join(names)

# update parameter dictionary 
opts["data_path"] = data_path

options = f"start time: {time.ctime()}\n" 
options += "  + command: python " + f"{' '.join(sys.argv)}\n"
options += "\n".join("  + {}: {}".format(k, v) for k, v in opts.items())

print(options)

if not args.background:
    print(f"do you wish to proceed?")
    proceed = input("[y/N]: ")

    if proceed.lower() == "n":
        print("you decided to not proceed. exiting.")
        exit()

pathlib.Path(f"{data_path}/logs").mkdir(parents=True, exist_ok=True)

# save run information to run log file, append if necessary
simInfoFile = f"{data_path}/logs/summary_simulation_info.log"
with open(simInfoFile, mode="a") as historyFile:
    historyFile.write(options)

print("submitting your jobs.")

if __name__=="__main__":
    mp.freeze_support()
    # do main process
    main(args.num_threads, collectParameters(opts), simInfoFile)
