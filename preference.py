import numpy as np
from scipy.integrate import odeint
from style import *
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import dash
from dash import dcc, html


def parse_params(cfus=False):
    if cfus:
        df = dict(pd.read_csv("parameters_cfus.csv"))
        params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
        p = params
        return p
    else:
        df = dict(pd.read_csv("parameters.csv"))
        params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
        p = params
        return p


def model(y, t, p):
    N1, N2, C1, C2 = y
    J1 =
