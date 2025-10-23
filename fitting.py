import pandas as pd
import numpy as np
import plotly.graph_objects as go
from os import path
from style import *
from scipy.integrate import odeint

d = "data/251018_succinate_glucose_plate_reader/metaod"


def parse_meta_od(dir):
    return pd.read_csv(path.join(dir, "metadata.csv")), pd.read_csv(
        path.join(dir, "measurements.csv"),
        index_col=0,
    )


class Condition:
    def __init__(self, species, carbon_source, cs_conc, signal, linegroups, data):
        self.species_names = species
        self.carbon_source = carbon_source
        self.concentration = cs_conc
        self.linegroups = linegroups
        self.signal = signal
        self.ys = []
        self.xs = []
        for lg in self.linegroups:
            x = data[lg + "_time"].to_numpy()
            y = data[lg + "_measurement"].to_numpy()
            self.xs.append(x)
            self.ys.append(y)
        self.y_average = np.mean(self.ys, axis=0)

    def filter_time(self, cut_off):
        filtered_xs = []
        filtered_ys = []
        for x, y in zip(self.xs, self.ys):
            mask = x <= cut_off
            filtered_xs.append(x[mask])
            filtered_ys.append(y[mask])
        self.xs = filtered_xs
        self.ys = filtered_ys

    def plot_condition(self):
        fname = f"{'+'.join(self.species_names)}_{self.carbon_source}_{self.concentration}mM_{self.signal}.svg"
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
        fig.update_layout(
            title=f"{'+'.join(self.species_names)} - {self.carbon_source} {self.concentration} mM",
            xaxis_title="Time (hours)",
            yaxis_title=self.signal,  # <- was "OD"
            width=600,
            height=400,
        )
        fig = style_plot(fig, font_size=14, line_thickness=2)
        return fig, fname


class Species:
    def __init__(self, name):
        self.name = name
        self.v_succ = None
        self.v_gluc = None
        self.K_succ = None
        self.K_gluc = None
        self.q_succ = None
        self.q_gluc = None
        self.t_lag = None
        self.a = None
        self.N0 = None
        self.y = None


import numpy as np
import pandas as pd
from scipy.integrate import odeint
import plotly.graph_objects as go

EPS = 1e-12


class Fitting:
    def __init__(self, condition):
        self.condition = condition
        self.species = []
        self.glucose = None
        self.succinate = None
        # average time grid (assumes aligned replicates)
        self.x = np.mean(self.condition.xs, axis=0)

        # --- species params ---
        df = pd.read_csv("parameters_species_fitting.csv", index_col=0)
        for species_name in condition.species_names:
            self.species.append(Species(species_name))

        for sp in self.species:
            p = df.loc[sp.name]
            sp.v_succ = float(p["v_succ"])
            sp.v_gluc = float(p["v_gluc"])
            sp.K_succ = float(p["K_succ"])
            sp.K_gluc = float(p["K_gluc"])
            sp.q_succ = float(p["q_succ"])
            sp.q_gluc = float(p["q_gluc"])
            sp.a = float(p["a"])
            sp.N0 = float(p["N0"])
            # NEW (Option B): acclimation timescale and initial activity
            sp.t_lag = float(p["lag"])
            sp.E0 = float(p["E0"]) if "E0" in p.index else 0.0

        # --- experiment / medium params ---
        df = pd.read_csv("parameters_experiment_fitting.csv")
        self.D = float(df["D"][0])
        self.M_succinate = float(df["succinate"][0])
        self.M_glucose = float(df["glucose"][0])

    # ---------- ORIGINAL (no lag) ----------
    def mac_arthur_mono_culture(self, y, t):
        sp = self.species[0]
        N, S, G = y

        # per-resource fluxes (no acclimation)
        J_S = sp.a * sp.v_succ * S / (sp.K_succ + S + EPS)
        J_G = (1.0 - sp.a) * sp.v_gluc * G / (sp.K_gluc + G + EPS)
        J = J_S + J_G

        dSdt = -(J_S * N) / sp.q_succ + self.D * (self.M_succinate - S)
        dGdt = -(J_G * N) / sp.q_gluc + self.D * (self.M_glucose - G)
        dNdt = (J * N) - self.D * N
        return [dNdt, dSdt, dGdt]

    def integrate_mac_arthur_mono_culture(self):
        y0 = [self.species[0].N0, self.M_succinate, self.M_glucose]
        N, S, G = odeint(self.mac_arthur_mono_culture, y0, self.x).T
        self.species[0].y = N
        self.succinate, self.glucose = S, G

        fig = self.condition.plot_condition()[0]
        fig.add_trace(
            go.Scatter(
                x=self.x,
                y=N,
                mode="lines",
                line=dict(width=2, color="red"),
                name="Fitted Model",
            )
        )
        return fig

    # ---------- OPTION B: acclimation state E ----------
    def mac_arthur_mono_culture_lag(self, y, t):
        sp = self.species[0]
        N, S, G = y

        FIXED_WIDTH = 1  # hours; change once here if needed
        k = 4.394 / FIXED_WIDTH  # 10–90% rise time mapping
        phi = 1.0 / (1.0 + np.exp(-k * (t - sp.t_lag)))

        EPS = 1e-12
        J_S = phi * sp.a * sp.v_succ * S / (sp.K_succ + S + EPS)
        J_G = phi * (1.0 - sp.a) * sp.v_gluc * G / (sp.K_gluc + G + EPS)
        J = J_S + J_G

        dSdt = -(J_S * N) / sp.q_succ + self.D * (self.M_succinate - S)
        dGdt = -(J_G * N) / sp.q_gluc + self.D * (self.M_glucose - G)
        dNdt = (J * N) - self.D * N
        return [dNdt, dSdt, dGdt]

    def integrate_mac_arthur_mono_culture_lag(self):
        sp = self.species[0]
        y0 = [sp.N0, self.M_succinate, self.M_glucose]
        N, S, G = odeint(self.mac_arthur_mono_culture_lag, y0, self.x).T
        sp.y = N
        self.succinate, self.glucose = S, G

        fig = self.condition.plot_condition()[0]
        fig.add_trace(
            go.Scatter(
                x=self.x,
                y=N,
                mode="lines",
                line=dict(width=2, color="red"),
                name="Fitted Model (lag)",
            )
        )
        return fig


def initialize_conditions(d):
    conditions = []
    meta, data = parse_meta_od(d)

    m = meta.copy()
    for col in ["exp_ID", "species", "carbon_source"]:
        m[col] = m[col].astype(str)

    required = ["exp_ID", "species", "carbon_source", "cs_conc", "linegroup"]
    m = m.dropna(subset=required)

    for (species, cs, conc), df in m.groupby(
        ["species", "carbon_source", "cs_conc"], dropna=False
    ):
        species_list = species.split("+") if "+" in species else [species]

        # Within each condition, build one Condition per experiment
        for experiment, g in df.groupby("exp_ID"):
            exp_l = str(experiment).lower()
            if exp_l.endswith("mcherry"):
                signal = "mcherry"
            elif exp_l.endswith("gfp"):
                signal = "gfp"
            else:
                signal = "OD"

            lgs = sorted(g["linegroup"].dropna().unique().tolist())
            conditions.append(Condition(species_list, cs, conc, signal, lgs, data))

    return conditions


_conditions_cache = None


def all_conditions(d):
    global _conditions_cache
    if _conditions_cache is None:
        _conditions_cache = initialize_conditions(d)
    return _conditions_cache


def get_condition(species, carbon_source, cs_conc, signal):
    for c in all_conditions(d):
        if (
            c.species_names == species
            and c.carbon_source == carbon_source
            and c.concentration == cs_conc
            and c.signal == signal
        ):
            return c


c = get_condition(["At"], "Succinate", 15, "OD")
c.filter_time(27.0)
F = Fitting(c)
fig = F.integrate_mac_arthur_mono_culture_lag()
# Glucose
# Oa	0.6	0.45	0.4	0.01	0.024	0.018	0			0.005	18
# At	0.28	1.1	0.5	0.01	0.06	0.058	0			0.007	10 with FIXED_WIDTH=10

fig.write_image("plots/fitting/plot.svg")

# Succinate
# At	0.32	1.1	0.5	0.01	0.03	0.058	1			0.007	2
