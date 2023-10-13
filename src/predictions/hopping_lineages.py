import numpy as np
import os
import matplotlib.pyplot as plt
from numpy import random
import argparse


def nearest_hotspot(hotspotlist, point, Lx, coeff=1, verbose=False):
    """
    find the closest hotspot downstream and within the parabola from (x,y)
    assumes hotspotlist is [[x1,y1],[x2,y2],...] and sorted in ys ascending
    """
    x, y = point

    def in_parabola(hs):
        delta = np.abs(hs[0] - x)
        if delta > Lx / 2:
            delta = Lx - delta
        # fix: added constant term in calculation from wolfram
        # return (coeff * delta ** 2 - (1/(4*coeff)))< (y - hs[1])
        return coeff * (delta) ** 2 < (y - hs[1])
        # return (coeff * delta) ** 2 < (y - hs[1])

    startidx = np.searchsorted(hotspotlist[:, 1], y) - 1

    if verbose:
        print(startidx, y, hotspotlist[startidx])
    if startidx < 0:
        return None
    while not in_parabola(hotspotlist[startidx]):
        if verbose:
            print("w", startidx, y, hotspotlist[startidx])
        startidx -= 1
        if startidx < 0:
            return None
    if verbose:
        print("y", startidx, y, hotspotlist[startidx])
    return hotspotlist[startidx]


def generate_hotspots(N, Lx, Ly):
    """poisson distributed points between [-Lx/2,Lx/2] and [0,Ly]"""
    xpts = random.random(N) * Lx
    ypts = random.random(N) * Ly
    return np.vstack((xpts, ypts)).T


def readFile(file):
    """read npz file, extract hotspots array"""
    # data = h5py.File(file, "r")
    data = np.load(file)
    # lineageheatmap = data["lineagesHeatmap"]
    hotspotslist = data["hsCarr"]
    x, y = zip(*[(point[0], point[1]) for point in hotspotslist])
    # return np.vstack((x, y)).T, lineageheatmap
    return np.vstack((x, y)).T


def trace_path(startingpoint, hotspotlist, Lx, coeff=1, root=False):
    """go downwards from startingpoint"""
    path = [startingpoint]
    # sort by y values
    hotspotlist = hotspotlist[hotspotlist[:, 1].argsort()]
    nextpoint = nearest_hotspot(hotspotlist, startingpoint, Lx, coeff)

    while nextpoint is not None:
        path.append(nextpoint)
        nextpoint = nearest_hotspot(hotspotlist, nextpoint, Lx, coeff)
    # insert a base point vertically below last point
    if root:
        path.append([path[-1][0], 0])
    return path


def plotParabola(
    pt, Lx, coeff=1, line=True, shade=False, n=50, lw=0.2, color="C0", alpha=0.1
):
    xpts = np.linspace(0, Lx, n)
    if line:
        plt.plot(xpts, pt[1] + coeff * (pt[0] - xpts) ** 2, "0.5", linewidth=lw)
    if shade:
        plt.fill_between(
            xpts,
            pt[1] + coeff * (pt[0] - xpts) ** 2,
            1000000000,
            color=color,
            alpha=alpha,
        )


def plot_path(points, Lx, fmt="yellow", lw=1):
    points = np.array(points)
    for i in range(len(points) - 1):
        if np.abs(points[i, 0] - points[i + 1, 0]) < Lx / 2:

            plt.plot(points[i : i + 2, 0], points[i : i + 2, 1], fmt, linewidth=lw)
        else:
            # plot two lines off sides, clip with limits later
            ycoords = points[i : i + 2, 1]
            xcoords = points[i : i + 2, 0]
            if xcoords[1] > xcoords[0]:
                plt.plot(xcoords + np.array([Lx, 0]), ycoords, fmt, linewidth=lw)
                plt.plot(xcoords - np.array([0, Lx]), ycoords, fmt, linewidth=lw)
            else:
                plt.plot(xcoords - np.array([Lx, 0]), ycoords, fmt, linewidth=lw)
                plt.plot(xcoords + np.array([0, Lx]), ycoords, fmt, linewidth=lw)


def createSim(coef, hotspots, sources, Lx, Ly, output, showParabola=False):
    fig, ax = plt.subplots(figsize=(12, 12))
    # plt.title(r"$y = %2.3f x^2$" % coeff)
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)

    ax.scatter(hotspots[:, 0], hotspots[:, 1], alpha=0.5)

    startends = dict()
    if showParabola:
        for pt in hotspots:
            plotParabola(pt, Lx, coef, shade=True, alpha=0.01, color="yellow")

    tmp_data = dict()
    trails = []
    for i in range(Lx):
        source = sources[i]
        backpath = trace_path(source, hotspots, Lx, coef, root=True)
        trails.append(backpath)
        plot_path(backpath, Lx)
        # plt.scatter([i], [Ly], marker="x", color="r")
        tmp_data[i] = backpath[-1]
    startends[coef] = tmp_data

    trails = np.array(trails, dtype=object)

    figpath = f"{output}/theory/parabolic/"
    if not os.path.isdir(figpath):
        os.makedirs(figpath)
    fig.savefig(f"{figpath}/hoppingLineages_{coef:.5f}.png")
    return startends


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="hopping lineages")
    parser.add_argument("--printInfo", action="store_true", help="enable printing")
    parser.add_argument("--output", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--coefficient", required=False, type=float)
    args = parser.parse_args()

    # task: implement readFile 
    lx, ly, intensity, density, radius, hs, sources = readFile(args.dataFile)

    if args.coefficient is None:
        estCoeff = (1 + intensity) / (2 * radius * intensity)
    else:
        estCoeff = args.coefficient

    print(f" we estimated coeff is {estCoeff}")

    createSim(estCoeff, hs, sources, lx, ly, args.output, showParabola=False)
