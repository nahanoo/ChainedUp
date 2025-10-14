from scipy.integrate import odeint
import numpy as np
import plotly.graph_objects as go
from style import *
import pandas as pd


meta = pd.read_csv("data/250918_glucose_lactose_plate/metadata.csv")
data = pd.read_csv("data/250918_glucose_lactose_plate/measurements.csv", index_col=0)


def get_linegroups(species, cs, signal, meta):
    linegroups = meta[
        (meta["species"] == species)
        & (meta["carbon_source"] == cs)
        & (meta["exp_ID"] == signal)
    ]["linegroup"].unique()
    return linegroups


def map_od_to_gfp(cs):
    species = "At"
    signal = "glucose_lactose_250918_glucose_lactose_screening_OD"
    lgs_od = get_linegroups(species, cs, signal, meta)
    signal = "glucose_lactose_250918_glucose_lactose_screening_GFP"
    lgs_gfp = get_linegroups(species, cs, signal, meta)
    ods = np.array([data[lg_od + "_measurement"].to_numpy() for lg_od in lgs_od]).mean(
        axis=0
    )[1:-1]
    gfps = np.array(
        [data[lg_gfp + "_measurement"].to_numpy() for lg_gfp in lgs_gfp]
    ).mean(axis=0)[1:-1]

    keep_od = []
    keep_gfp = []

    last_od, last_gfp = -np.inf, -np.inf

    for od, gfp in zip(ods, gfps):
        if (od > last_od) and (gfp > last_gfp):
            keep_od.append(od)
            keep_gfp.append(gfp)
            last_od, last_gfp = od, gfp

    gfp_to_od = pd.DataFrame({"od": keep_od, "gfp": keep_gfp})
    gfp_to_od.to_csv("data/gfp_to_od.csv", index=False)

    species = "Oa"
    signal = "glucose_lactose_250918_glucose_lactose_screening_OD"
    lgs_od = get_linegroups(species, cs, signal, meta)
    signal = "glucose_lactose_250918_glucose_lactose_screening_mcherry"
    lgs_mcherry = get_linegroups(species, cs, signal, meta)
    ods = np.array([data[lg_od + "_measurement"].to_numpy() for lg_od in lgs_od]).mean(
        axis=0
    )[1:-1]
    gmcherrys = np.array(
        [data[lg_mcherry + "_measurement"].to_numpy() for lg_mcherry in lgs_mcherry]
    ).mean(axis=0)[1:-1]
    keep_od = []
    keep_mcherry = []

    last_od, last_mcherry = -np.inf, -np.inf

    for od, mcherry in zip(ods, gmcherrys):
        if (od > last_od) and (mcherry > last_mcherry):
            keep_od.append(od)
            keep_mcherry.append(mcherry)
            last_od, last_mcherry = od, mcherry

    gfp_to_od = pd.DataFrame({"od": keep_od, "mcherry": keep_mcherry})
    gfp_to_od.to_csv("data/mcherry_to_od.csv", index=False)
