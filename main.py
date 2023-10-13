import argparse
from numpy import floor, arange, array
import pathlib
import os
import time
import multiprocessing as mp
import sys
import copy
from random import random, seed, randint


class KeyboardInterruptError(Exception):
    pass


class Parameters:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def collectParameters(cmdinput) ->list[Parameters]:
    parameterList = []
    start, step, stop = cmdinput["intervals"]
    
    # initialize rng with cli seed
    seed(cmdinput["rngSeed"])

    # create multiple collections for each landscape config
    for ithC in range(1, cmdinput["nCollections"] + 1):
        # ! for each collection, set a random seed
        rng = randint(1, 1e10)

        # step through range of parameter parameter values 
        for parameterVal in arange(start, stop + step, step=step):
            cdict = copy.deepcopy(cmdinput)

            # modify scanned parameter and savepath
            cdict[cmdinput["parameter"]] = parameterVal
            parameters = Parameters(**cdict)
            parameters.savepath = f'{cdict["savepath"]}/collection-{ithC}'
            parameters.rngSeed = rng

            parameterList.append((parameters))
    return parameterList


def to_list(arg):
    return [float(i) for i in arg.split(",")]


def runJulia(prms: Parameters):
    cmd_roughfront = "--roughfront" if prms.roughfront else ""
    cmd_overwrite = "--overwrite" if prms.overwrite else ""
    cmd_landscape = f"--landscape {prms.landscape}" if prms.landscape is not None else ""
    cmd_hsRegion = f"--hsRegion {int(prms.hsRegion[0])},{int(prms.hsRegion[1])}" if prms.hsRegion is not None else ""

    command = f"""julia --project=@. -O3 main.jl \
            --numberTrials {prms.numberTrials} \
            --numberSamples {int(prms.numberSamples)} \
            --dims {prms.dims[0]},{prms.dims[1]} \
            --density {prms.density} \
            --radius {prms.radius} \
            --intensity {prms.intensity} \
            {cmd_hsRegion} \
            --rngSeed {int(prms.rngSeed)} \
            --savePath {prms.savepath} \
            --parameter {prms.parameter} \
            --sampleImages \
            {cmd_landscape} \
            {cmd_overwrite} \
            {cmd_roughfront} \
            1> /dev/null
        """
    print(f"starting: (D={round(prms.density,3)}, R={prms.radius}, I={prms.intensity})")
    os.system(command)
    print(f"\nfinished: (D={round(prms.density,3)}, R={prms.radius}, I={prms.intensity})")


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

    msg = f"finished. Elapsed wall time: {round((time.time()-t1)/3600, 3)} hour(s), at {time.ctime()}"
    print(msg)

    with open(simInfoFile, mode="a") as historyFile:
        historyFile.write(msg)

parser = argparse.ArgumentParser()
parser.add_argument("--numberTrials", type=int, default=1)
parser.add_argument("--numberSamples", type=int, default=1000)
parser.add_argument("--dims", required=True, type=to_list)
parser.add_argument("--density", type=float, default=0.0)
parser.add_argument("--radius", type=int, default=10)
parser.add_argument("--intensity", type=float, default=3.0)
parser.add_argument("--hsRegion", type=to_list)
parser.add_argument("--rngSeed", type=int, default=1)
parser.add_argument("--overwrite", action="store_true")
parser.add_argument("--savepath", required=True)
parser.add_argument("--landscape", default=None)

parser.add_argument("--parameter", default="density")
parser.add_argument("--intervals", required=True, type=to_list)
parser.add_argument("--nCollections", type=int, default=1)
parser.add_argument("--num_threads", type=int, default=2)
parser.add_argument("--background", action="store_true")
parser.add_argument("--roughfront", action="store_true")

args = parser.parse_args()
cmdinput = vars(args)

cmdinput[args.parameter] = "xxx"
pc = Parameters(**cmdinput)

names = [
    f"dims:{pc.dims[0]},{pc.dims[1]}",
    f'hsRegion:{pc.hsRegion}',
    f'I:{pc.intensity}',
    f'D:{pc.density}',
    f'R:{pc.radius}',
    f"trials:{pc.numberTrials}",
    f"nCollections:{pc.nCollections}",
    f'intervals:{",".join([str(s) for s in pc.intervals])}',
    f'landscape:{"fixed" if pc.landscape is not None else "varied"}']

savepath = f"{args.savepath}/" + "_".join(names)

# update parameter dictionary 
cmdinput["savepath"] = savepath

options = f"start time: {time.ctime()}\n" 
options += "  + command: python " + f"{' '.join(sys.argv)}\n"
options += "\n".join("  + {}: {}".format(k, v) for k, v in cmdinput.items())

print(options)

if not args.background:
    print(f"do you wish to proceed?")
    proceed = input("[y/N]: ")

    if proceed.lower() == "n":
        print("you decided to not proceed. exiting.")
        exit()

pathlib.Path(f"{savepath}/logs").mkdir(parents=True, exist_ok=True)

# save run information to run log file, append if necessary
simInfoFile = f"{savepath}/logs/summary_simulation_info.log"
with open(simInfoFile, mode="a") as historyFile:
    historyFile.write(options)

print("submitting your jobs.")

# for c in collectParameters(cmdinput):
#     print("\n")
#     for k, v in vars(c).items():
#         print(k, " : ", v)

# exit()

if __name__=="__main__":
    mp.freeze_support()
    # do main process
    main(args.num_threads, collectParameters(cmdinput), simInfoFile)
