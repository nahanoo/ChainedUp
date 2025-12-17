from statistics import stdev, mean
from os.path import join
import pandas as pd
from glob import glob
import json
import numpy as np
from sympy import *


def cfu_parser(e):
    with open(join("/", "home", "eric", "ChainedUp", "data", e, "order.txt"), "r") as f:
        order = f.read().rstrip().split(",")
    """Parses counted CFUs based on template xlsx"""
    # For eveery species there is one xlsx sheet
    fs = glob(join("/", "home", "eric", "ChainedUp", "data", e, "cfus*.xlsx"))
    # df for concatenating species sheets
    dfs = []
    for f in fs:
        # Species name
        species = f[-7:-5]
        cfus = pd.read_excel(f, header=1)
        cfus.insert(0, "species", species)
        dfs.append(cfus)
    # Reindexing because of concat
    df = pd.concat(dfs)
    df.index = range(len(df))

    # Splitting "|" separated triplicates
    counts = []
    for count, dilution in zip(df["count"], df["dilution"]):
        # Conversion to CFUs/mL sample volume 5 uL
        counts.append([(int(n) / 10**dilution) * 200 for n in count.split("|")])
    # Calculating average and stdev
    df["average"] = [mean(count) for count in counts]
    df["stdev"] = [stdev(count) for count in counts]

    # Adding total counts for composition
    df.insert(len(df.columns), "total", None)
    for t in df["sample_time"]:
        # Subsetting for reactor and sample time
        for reactor in order:
            tmp = df[df["sample_time"] == t]
            tmp = tmp[tmp["reactor"] == reactor]
            total = tmp["average"].sum()
            for i in tmp.index:
                df.at[i, "total"] = total
    # Adding composition
    df.insert(len(df.columns), "composition", None)
    df["composition"] = df["average"] / df["total"]
    return df, order


def chibio_parser(e, down_sample=False, sample_size=10):

    def average(df):
        out = pd.DataFrame(
            columns=["exp_time", "od_measured", "FP1_base", "FP1_emit1", "FP1_emit2"]
        )
        i = 0
        j = 0
        while True:
            try:
                sliece = df.loc[i : i + sample_size]
                out.loc[j] = [
                    sliece.iloc[-1]["exp_time"],
                    np.average(sliece["od_measured"]),
                    np.average(sliece["FP1_base"]),
                    np.average(sliece["FP1_emit1"]),
                    np.average(sliece["FP1_emit2"]),
                ]
                i += sample_size
                j += 1
            except IndexError:
                break
        return out

    with open(join("/", "home", "eric", "ChainedUp", "data", e, "order.txt"), "r") as f:
        order = f.read().rstrip().split(",")
    dfs = []
    for reactor in order:
        f = glob(
            join("/", "home", "eric", "ChainedUp", "data", e, reactor, "*data.csv")
        )[0]
        data = pd.read_csv(
            f, usecols=["exp_time", "od_measured", "FP1_base", "FP1_emit1", "FP1_emit2"]
        )
        if down_sample:
            data = average(data)
        data.insert(1, "reactor", reactor)
        data["exp_time"] = data["exp_time"] / 60 / 60
        dfs.append(data)

    df = pd.concat(dfs)
    dfs = []
    for c in ["od_measured", "FP1_emit1", "FP1_emit2", "FP1_base"]:
        tmp = pd.DataFrame(columns=["exp_time", "reactor", "measurement", "sensor"])
        tmp["reactor"] = df["reactor"]
        tmp["exp_time"] = df["exp_time"]
        tmp["measurement"] = df[c]
        tmp["sensor"] = c
        dfs.append(tmp)
    df = pd.concat(dfs)

    return df, order


def fluorescence_paresr(e):
    with open(join("/", "home", "eric", "ChainedUp", "data", e, "order.txt"), "r") as f:
        order = f.read().rstrip().split(",")
    dfs = []
    for reactor in order:
        f = glob(
            join("/", "home", "eric", "ChainedUp", "data", e, reactor, "*data.csv")
        )[0]
        columns = [
            "exp_time",
            "od_measured",
            "FP1_base",
            "FP1_emit1",
            "FP1_emit2",
            "FP2_base",
            "FP2_emit1",
            "FP2_emit2",
            "FP3_base",
            "FP3_emit1",
            "FP3_emit2",
            "reactor",
            "raw_FP1_emit1",
            "raw_FP1_emit2",
            "raw_FP2_emit1",
            "raw_FP2_emit2",
            "raw_FP3_emit1",
            "raw_FP3_emit2",
        ]
        data_columns = columns[:11]
        data = pd.read_csv(
            f,
            usecols=[
                "exp_time",
                "od_measured",
                "FP1_base",
                "FP1_emit1",
                "FP1_emit2",
                "FP2_base",
                "FP2_emit1",
                "FP2_emit2",
                "FP3_base",
                "FP3_emit1",
                "FP3_emit2",
            ],
        )
        df = pd.DataFrame(columns=columns)
        for c in data_columns:
            df[c] = data[c]
        FP1_emit = columns[3:5]
        FP1_emit_raw = columns[12:14]
        for emit, raw in zip(FP1_emit, FP1_emit_raw):
            df[raw] = data[emit] * data["FP1_base"]

        FP2_emit = columns[6:8]
        FP2_emit_raw = columns[14:16]
        for emit, raw in zip(FP2_emit, FP2_emit_raw):
            df[raw] = data[emit] * data["FP2_base"]

        FP3_emit = columns[9:11]
        FP3_emit_raw = columns[16:18]
        for emit, raw in zip(FP3_emit, FP3_emit_raw):
            df[raw] = data[emit] * data["FP3_base"]

        df["reactor"] = reactor
        df["exp_time"] = df["exp_time"] / 60 / 60
        OD = df["od_measured"].to_numpy()
        window_size = 10
        OD = np.convolve(OD, np.ones(window_size), "same") / window_size
        df.at[:, "od_measured"] = OD
        dfs.append(df)

    out = pd.concat(dfs)

    return out


def ct_oa_parser(e, parse_cfus=True):
    df = fluorescence_paresr(e)
    Ms = [df[df["reactor"] == M] for M in sorted(set(df["reactor"]))]
    columns = [
        "exp_time",
        "reactor",
        "od_measured",
        "raw_FP1_emit1",
        "raw_FP2_emit1",
        "raw_FP2_emit2",
        "FP3_emit1",
    ]
    Ms = [M[columns] for M in Ms]
    od_history = 60
    rolling_window = 15
    for M in Ms:
        ods = []
        for i, row in M.iterrows():
            od = row["od_measured"]
            if i - od_history > 0:
                avg = np.average(M.loc[i - od_history : i, "od_measured"])
            else:
                avg = np.average(M.loc[:i, "od_measured"])
            if 1 / avg * od > 0.95:
                ods.append(od)
            else:
                ods.append(avg)
        OD = df["od_measured"].to_numpy()
        window_size = 10
        OD = np.convolve(OD, np.ones(window_size), "same") / window_size
        M.at[:, "od_measured"] = OD
    Ms = [M.iloc[::10] for M in Ms]
    if parse_cfus:
        cfus = cfu_parser(e)[0]
        cfus = {
            "ct": [
                cfus[(cfus["reactor"] == M) & (cfus["species"] == "ct")]
                for M in sorted(set(df["reactor"]))
            ],
            "oa": [
                cfus[(cfus["reactor"] == M) & (cfus["species"] == "oa")]
                for M in sorted(set(df["reactor"]))
            ],
        }
    else:
        cfus = None
    return (Ms, cfus)


def calibration():
    # This function was used before calibration_csv so it's not used anymore
    # but remains for old data
    # Known data points for OD and R
    R1, R2, OD1, OD2, k, b = symbols("R1 R2 OD1 OD2 k b")

    # Set up the adjusted equation with m fixed
    eq1 = Eq(OD1, k * log(R1, 10) + b)
    eq2 = Eq(OD2, k * log(R2, 10) + b)

    # Solve for k and b
    solution = solve((eq1, eq2), (k, b))
    return solution[k], solution[b]


def calibration_csv(csv, data_fluorescence_parser):
    R0, R1, OD0, OD1, k, b = symbols("R1 R2 OD1 OD2 k b")

    # Set up the adjusted equation with m fixed
    eq1 = Eq(OD0, k * log(R0, 10) + b)
    eq2 = Eq(OD1, k * log(R1, 10) + b)

    # Solve for k and b
    solution = solve((eq1, eq2), (k, b))
    k_eq, b_eq = solution[k], solution[b]

    cal = pd.read_csv(csv, index_col="reactor")
    data = data_fluorescence_parser
    out = []
    for r in set(data["reactor"]):
        df = data[data["reactor"] == r]
        k = float(
            k_eq.subs(
                {
                    OD0: cal.loc[r]["OD0"],
                    OD1: cal.loc[r]["OD1"],
                    R0: cal.loc[r]["R0"],
                    R1: cal.loc[r]["R1"],
                }
            )
        )
        b = b_eq.subs(
            {
                OD0: cal.loc[r]["OD0"],
                OD1: cal.loc[r]["OD1"],
                R0: cal.loc[r]["R0"],
                R1: cal.loc[r]["R1"],
            }
        )
        OD = k * np.log10(df["od_measured"]) + b
        df.insert(len(df.columns), "od_calibrated", np.float64(OD))
        out.append(df)
    return pd.concat(out)
