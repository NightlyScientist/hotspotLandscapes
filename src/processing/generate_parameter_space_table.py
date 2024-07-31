# %%
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True)
parser.add_argument("--processed", action="store_true")
parser.add_argument("--deprecated", action="store_true")
args = parser.parse_args()
# opts = vars(args)

basePath = args.input

# %%
df = None

for path in os.listdir(basePath):
    if not ("XY:" in path or "dims" in path):
        continue

    for subpath in os.listdir(os.path.join(basePath, path)):
        if "logs" in subpath:
            continue

        for sspath in os.listdir(os.path.join(basePath, path, subpath)):
            _opts = os.path.join(basePath, path, subpath, sspath, "Opts.csv")
            opts = pd.read_csv(_opts, header=0, delimiter=";").squeeze("rows").to_dict()
            opts["path"] = os.path.join(basePath, path, subpath, sspath)

            dct = {k: [v] for k, v in opts.items()}

            if df is None:
                df = pd.DataFrame.from_dict(dct)
            else:
                _df = pd.DataFrame.from_dict(dct)
                df = pd.concat([df, _df], ignore_index=True)

# %%
# group entries
groupParameters = ["ref_line", "height", "width", "radius", "gap"]
abbr = ["rf", "H", "W", "R", "G"]
gdf = df.groupby(groupParameters)


# %%
def partition(x, threshold=0.002):
    _x = np.sort(np.array(x))
    _indx = np.argsort(np.array(x))
    cuts = np.where(np.diff(_x) > threshold)
    x_p = np.split(_x, cuts[0] + 1)
    i_p = np.split(_indx, cuts[0] + 1)
    return x_p, i_p


def processed(df):
    dct = {}
    for i in range(len(df)):
        p = df.loc[i, "path"]
        density = df.loc[i, "density"]
        intensity = df.loc[i, "intensity"]
        for k in filter(lambda k: k.startswith("_"), os.listdir(p)):
            key = k.split(".")[0]
            t = (density, intensity)
            if not key in dct:
                dct[key] = [t]
            dct[key].append(t)
    return dct


# %%
flag = "_deprecated" if args.deprecated else ""
sb.set_theme("paper")

# .save parameter table in the same folder as the collection folder
savePath = args.input + "/logs/"

for key in gdf.groups.keys():
    ps = gdf.get_group(key).reset_index(drop=True)

    _params = "_".join([f"{k}:{v}" for (k, v) in zip(abbr, key)])

    indx_groups = partition(ps.density)[1]
    new_df = ps.copy()

    if not args.processed:
        new_df["group"] = np.zeros(len(new_df.density), dtype=int)
        for g in range(len(indx_groups)):
            for i in indx_groups[g]:
                new_df.loc[i, "group"] = g

        # save new dataframe
        new_df.to_csv(savePath + f"/PS_{_params}{flag}.csv", sep=";")

        # save fig of parameter space and grouped dataframe
        fig, ax = plt.subplots()
        fig.suptitle(f"Sampled Parameter Space\n{_params}")
        sb.scatterplot(data=ps, x="density", y="intensity", marker="o", ax=ax)
        fig.savefig(savePath + f"/PSView_{_params}{flag}.png")

    else:
        # save fig of process metrics
        for k, xy in processed(ps).items():
            fig, ax = plt.subplots()
            fig.suptitle(f"Sampled {k} in Parameter Space")
            sb.scatterplot(data=pd.DataFrame(xy), x=0, y=1, ax=ax)
            fig.savefig(savePath + f"/{k}_PSView_{_params}{flag}.png")

# %%
