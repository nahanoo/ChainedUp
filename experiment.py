import numpy as np
import pandas as pd
from os import path
import plotly.graph_objects as go
from scipy.integrate import odeint
from style import style_plot


class Species:
    def __init__(self, name: str, row: pd.Series):
        self.name = name
        self.v_succ = float(row["v_succ"])
        self.v_gluc = float(row["v_gluc"])
        self.K_succ = float(row["K_succ"])
        self.K_gluc = float(row["K_gluc"])
        self.q_succ = float(row["q_succ"])
        self.q_gluc = float(row["q_gluc"])
        self.t_lag = float(row["t_lag"])
        self.a = float(row["a"])
        self.N0 = float(row["N0"])
        self.lag = self.t_lag
        self.y = None


class Condition:
    def __init__(
        self, species, carbon_source, cs_conc, signal, linegroups, data, species_df
    ):
        self.species_names = list(species)
        self.species = []
        self.carbon_source = carbon_source
        self.concentration = cs_conc
        self.linegroups = list(linegroups)
        self.signal = signal

        self.xs = [data[f"{lg}_time"].to_numpy() for lg in self.linegroups]
        self.ys = [data[f"{lg}_measurement"].to_numpy() for lg in self.linegroups]

        for s in self.species_names:
            row = species_df.loc[s]
            self.species.append(Species(s, row))

    def filter_time(self, cut_off: float):
        new_xs, new_ys = [], []
        for x, y in zip(self.xs, self.ys):
            m = x <= cut_off
            new_xs.append(x[m])
            new_ys.append(y[m])
        self.xs, self.ys = new_xs, new_ys

    def plot_condition(self):

        fig = go.Figure()
        for x, y in zip(self.xs, self.ys):
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    line=dict(width=1, color="purple"),
                    showlegend=False,
                )
            )
        ylabels = {"OD": "OD600", "gfp": "GFP (a.u.)", "mcherry": "mCherry (a.u.)"}
        fig.update_layout(
            title=f"{'+'.join(self.species_names)} - {self.carbon_source} {self.concentration} mM",
            xaxis_title="Time (hours)",
            yaxis_title=ylabels.get(self.signal, self.signal),
            width=600,
            height=400,
        )
        fig = style_plot(fig, font_size=14, line_thickness=2)
        return fig


class Experiment:
    def __init__(
        self,
        d,
        species_param_path="parameters/parameters_species_fitting.csv",
        exp_param_path="parameters/parameters_experiment_fitting.csv",
    ):
        self.meta = pd.read_csv(path.join(d, "metadata.csv"))
        self.data = pd.read_csv(path.join(d, "measurements.csv"), index_col=0)
        self.species_df = pd.read_csv(species_param_path, index_col=0)

        exp_df = pd.read_csv(exp_param_path)
        self.D = float(exp_df["D"][0])
        self.M_succinate = float(exp_df["succinate"][0])
        self.M_glucose = float(exp_df["glucose"][0])

        self.conditions = []
        self.build_conditions()

    def build_conditions(self):
        m = self.meta.copy()
        for col in ["exp_ID", "species", "carbon_source"]:
            m[col] = m[col].astype(str)
        required = ["exp_ID", "species", "carbon_source", "cs_conc", "linegroup"]
        m = m.dropna(subset=required)

        for (species, cs, conc), df in m.groupby(
            ["species", "carbon_source", "cs_conc"], dropna=False
        ):
            species_list = species.split("+") if "+" in species else [species]
            for experiment, g in df.groupby("exp_ID"):
                exp_l = str(experiment).lower()
                if exp_l.endswith("mcherry"):
                    signal = "mcherry"
                elif exp_l.endswith("gfp"):
                    signal = "gfp"
                else:
                    signal = "OD"

                lgs = sorted(g["linegroup"].dropna().unique().tolist())
                self.conditions.append(
                    Condition(
                        species_list, cs, conc, signal, lgs, self.data, self.species_df
                    )
                )

    def get_condition(self, species, carbon_source, cs_conc, signal):
        target = (tuple(sorted(species)), carbon_source, cs_conc, signal)
        for c in self.conditions:
            key = (
                tuple(sorted(c.species_names)),
                c.carbon_source,
                c.concentration,
                c.signal,
            )
            if key == target:
                return c
        return None
