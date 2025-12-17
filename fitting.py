from experiment import Species, Experiment
from os import path
import pandas as pd
import numpy as np
from style import *
from scipy.integrate import odeint
import plotly.graph_objects as go


class Model:
    def __init__(self, at, oa, e, t, C_mono, D):
        self.at = at
        self.oa = oa
        self.e = e
        self.t = t
        self.C_mono = C_mono
        self.D = D

        self.C_to_mM_succinate = {
            45: 11.25,
            30: 7.5,
            15: 3.75,
            7.5: 1.875,
            5: 1.25,
            2.5: 0.625,
        }
        self.C_to_mM_glucose = {
            45: 7.5,
            30: 5,
            15: 2.5,
            7.5: 1.25,
            5: 0.833,
            2.5: 0.417,
        }
        self.C_to_mM_glucose_succinate = {
            90: [11.25, 7.5],
            60: [7.5, 5],
            30: [3.75, 2.5],
            15: [1.875, 1.25],
            10: [1.25, 0.833],
            5: [0.625, 0.417],
        }
        self.C_mono = C_mono
        self.M_succinate = self.C_to_mM_succinate[self.C_mono]
        self.M_glucose = self.C_to_mM_glucose[self.C_mono]
        self.glucose = None
        self.succinate = None
        self.D = D

    def lag_time(self, t, t_lag, lag_width, *, min_width=1e-6, clip=50.0, frac=0.95):

        # Handle trivial "no lag" case; keep broadcasting behavior
        if t_lag <= 0:
            return np.ones_like(t, dtype=float) if np.ndim(t) else 1.0

        # Left width and (steep) right width
        wl = max(float(lag_width), float(min_width))
        right_frac = float(getattr(self, "lag_right_frac", frac))  # tune if needed
        wr = max(wl * max(right_frac, 0.0), float(min_width))

        C_10_90 = 4.39444872453601  # 2*ln(9)
        k_left = C_10_90 / wl
        k_right = C_10_90 / wr

        x = np.asarray(t, dtype=float) - float(t_lag)

        # piecewise logistic with separate slopes, clipped for stability
        phi = np.empty_like(x, dtype=float)
        mask_left = x < 0
        if np.any(mask_left):
            xL = np.clip(k_left * x[mask_left], -clip, 0.0)
            phi[mask_left] = 1.0 / (1.0 + np.exp(-xL))
        if np.any(~mask_left):
            xR = np.clip(k_right * x[~mask_left], 0.0, clip)
            phi[~mask_left] = 1.0 / (1.0 + np.exp(-xR))

        # Return scalar if scalar input
        return phi if np.ndim(t) else float(phi)

    def mac_arthur_lag(self, y, t):
        EPS = 1e-12
        at, oa, g, s = y

        # Enforce non-negativity in the kinetics
        g_eff = 0.0 if g <= 0.0 else g
        s_eff = 0.0 if s <= 0.0 else s

        # ON/OFF gates (as you had)
        phi_S_at = (
            1
            if self.at.lag_succ == 0
            else self.lag_time(t, self.at.lag_succ, self.at.lag_succ_width, frac=0.95)
        )
        phi_G_at = (
            1
            if self.at.lag_gluc == 0
            else self.lag_time(t, self.at.lag_gluc, self.at.lag_gluc_width, frac=0.05)
        )
        phi_G_oa = (
            1
            if self.oa.lag_gluc == 0
            else self.lag_time(t, self.oa.lag_gluc, self.oa.lag_gluc_width)
        )
        phi_S_oa = (
            1
            if self.oa.lag_succ == 0
            else self.lag_time(t, self.oa.lag_succ, self.oa.lag_succ_width)
        )

        # Uptake fluxes use clamped resources
        J_S_at = (
            phi_S_at
            * self.at.a
            * self.at.v_succ_lag
            * s_eff
            / (self.at.K_succ + s_eff + EPS)
        )
        J_G_at = (
            phi_G_at
            * (1.0 - self.at.a)
            * self.at.v_gluc_lag
            * g_eff
            / (self.at.K_gluc + g_eff + EPS)
        )
        J_S_oa = (
            phi_S_oa
            * self.oa.a
            * self.oa.v_succ_lag
            * s_eff
            / (self.oa.K_succ + s_eff + EPS)
        )
        J_G_oa = (
            phi_G_oa
            * (1.0 - self.oa.a)
            * self.oa.v_gluc_lag
            * g_eff
            / (self.oa.K_gluc + g_eff + EPS)
        )

        J_at = J_S_at + J_G_at
        J_oa = J_S_oa + J_G_oa

        # D may be zero in batch; keep form general
        datdt = at * (J_at - self.D)
        doadt = oa * (J_oa - self.D)

        dgdt = self.D * (self.M_glucose - g_eff) - (
            J_G_at * at / self.at.q_gluc + J_G_oa * oa / self.oa.q_gluc
        )
        dsdt = self.D * (self.M_succinate - s_eff) - (
            J_S_at * at / self.at.q_succ + J_S_oa * oa / self.oa.q_succ
        )

        return [datdt, doadt, dgdt, dsdt]

    def calculate_lag(self):
        t = self.t
        at_phi_succ = np.array(
            [
                self.lag_time(ti, self.at.lag_succ, self.at.lag_succ_width, frac=0.05)
                for ti in t
            ]
        )
        at_phi_gluc = np.array(
            [
                self.lag_time(ti, self.at.lag_gluc, self.at.lag_gluc_width, frac=0.05)
                for ti in t
            ]
        )
        oa_phi_succ = np.array(
            [self.lag_time(ti, self.oa.lag_succ, self.oa.lag_succ_width) for ti in t]
        )
        oa_phi_gluc = np.array(
            [self.lag_time(ti, self.oa.lag_gluc, self.oa.lag_gluc_width) for ti in t]
        )
        self.at.phi_succ = at_phi_succ
        self.at.phi_gluc = at_phi_gluc
        self.oa.phi_succ = oa_phi_succ
        self.oa.phi_gluc = oa_phi_gluc

    def get_parameters(self):

        phi_S_at = (
            1
            if self.at.lag_succ == 0
            else self.lag_time(self.t, self.at.lag_succ, self.at.lag_succ_width)
        )
        phi_G_at = (
            1
            if self.at.lag_gluc == 0
            else self.lag_time(self.t, self.at.lag_gluc, self.at.lag_gluc_width)
        )
        phi_G_oa = (
            1
            if self.oa.lag_gluc == 0
            else self.lag_time(self.t, self.oa.lag_gluc, self.oa.lag_gluc_width)
        )
        phi_S_oa = (
            1
            if self.oa.lag_succ == 0
            else self.lag_time(self.t, self.oa.lag_succ, self.oa.lag_succ_width)
        )
        J_S_at = (
            phi_S_at
            * self.at.a
            * self.at.v_succ_lag
            * self.succinate
            / (self.at.K_succ + self.succinate)
        )
        J_G_at = (
            phi_G_at
            * (1.0 - self.at.a)
            * self.at.v_gluc_lag
            * self.glucose
            / (self.at.K_gluc + self.glucose)
        )

        J_S_oa = (
            phi_S_oa
            * self.oa.a
            * self.oa.v_succ_lag
            * self.succinate
            / (self.oa.K_succ + self.succinate)
        )
        J_G_oa = (
            phi_G_oa
            * (1.0 - self.oa.a)
            * self.oa.v_gluc_lag
            * self.glucose
            / (self.oa.K_gluc + self.glucose)
        )
        self.at.v_succ = np.max(J_S_at)
        self.at.v_gluc = np.max(J_G_at)
        self.oa.v_succ = np.max(J_S_oa)
        self.oa.v_gluc = np.max(J_G_oa)

    def integrate_mac_arthur(self):
        y0 = self.at.N0, self.oa.N0, self.M_glucose, self.M_succinate
        self.at.y, self.oa.y, self.glucose, self.succinate = odeint(
            self.mac_arthur_lag, y0, self.t
        ).T

    def plot_at_oa(self, fig=False):
        if not fig:
            fig = go.Figure()
        species = []
        if self.at.y[0] != 0:
            fig.add_trace(
                go.Scatter(
                    x=self.t,
                    y=self.at.y,
                    mode="lines",
                    line=dict(width=2, color="green"),
                    name="Model At",
                )
            )
            species.append("At")
        if self.oa.y[0] != 0:
            fig.add_trace(
                go.Scatter(
                    x=self.t,
                    y=self.oa.y,
                    mode="lines",
                    line=dict(width=2, color="red"),
                    name="Model Oa",
                )
            )
            species.append("Oa")
        if self.oa.y[0] != 0 and self.at.y[0] != 0:
            fig.add_trace(
                go.Scatter(
                    x=self.t,
                    y=self.at.y + self.oa.y,
                    mode="lines",
                    line=dict(width=2, color="purple"),
                    name="Model Total",
                )
            )
        title = f"{' + '.join(species)} - {self.M_glucose} mM Glucose, {self.M_succinate} mM Succinate"
        if not fig:
            fig.update_layout(
                title=title,
                xaxis_title="Time (hours)",
                yaxis_title="OD600",
                width=600,
                height=400,
            )
            fig = style_plot(fig, font_size=14, line_thickness=2)
        return fig


def simulate_growth_curves():
    e = Experiment(d="data/251018_succinate_glucose_plate_reader/metaod/")
    concentrations = [2.5, 5, 7.5, 15]
    carbon_sources = ["Succinate", "Glucose"]
    for conc in concentrations:
        # load parameter file once per concentration
        p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
        params = pd.read_csv(p_f, index_col=0)

        for cs in carbon_sources:
            a_pref = 1.0 if cs == "Succinate" else 0.0

            for focal in ("Oa", "At"):
                # build fresh species objects from table
                at = Species("At", params.loc["At"])
                oa = Species("Oa", params.loc["Oa"])

                # monoculture: zero the non-focal initial biomass
                if focal == "At":
                    oa.N0 = 0.0
                    at.a = a_pref
                else:
                    at.N0 = 0.0
                    oa.a = a_pref

                cond = e.get_condition([focal], cs, conc, "OD")
                if cond is None:
                    print(f"[WARN] Condition not found: {focal}, {cs}, {conc} mM, OD")
                    continue

                cond.filter_time(30)
                fig = cond.plot_condition()
                # integrate model on this condition's time grid
                tgrid = np.linspace(0, cond.xs[0][-1], 10000)
                M = Model(at, oa, e, tgrid, conc, 0)
                M.integrate_mac_arthur()
                M.plot_at_oa(fig)
                M.calculate_lag()
                M.get_parameters()
                sp = at if focal == "At" else oa
                lag_period = sp.phi_gluc if cs == "Glucose" else sp.phi_succ

                for i, p in enumerate(lag_period):
                    if p >= 0.95:
                        break
                fig.add_vrect(
                    x0=0,
                    x1=tgrid[i],
                    fillcolor="lightgray",
                    opacity=0.5,
                    line_width=0,
                    annotation_text="Lag period",  # your label
                    annotation_position="top",  # ⟵ centered inside the band
                )

                # print(sp.v_gluc)
                if cs == "Succinate":
                    mu_val, K_val = sp.v_succ, sp.K_succ
                else:
                    mu_val, K_val = sp.v_gluc, sp.K_gluc
                C_to_mM = M.C_to_mM_glucose if cs == "Glucose" else M.C_to_mM_succinate
                title = f"{focal} on {cs} {C_to_mM[conc]:g} mM — μ: {mu_val:g} 1/h   K: {K_val:g} mM"
                fig.update_layout(title=title)
                fig = style_plot(fig, font_size=14, line_thickness=2)

                safe_conc = f"{conc:g}".replace(".", "_")
                outfile = path.join(
                    "plots",
                    "fitting",
                    f"{focal.lower()}_{safe_conc}_mM_{cs.lower()}.svg",
                )
                fig.write_image(outfile)


def fit_at_glucose_succinate():
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)
    e = Experiment(d="data/251018_succinate_glucose_plate_reader/metaod/")
    e.build_conditions()
    xs = np.linspace(0, 28, 1000)
    succ_at = e.get_condition(["At"], "Succinate", conc, "OD")
    succ_at.filter_time(xs[-1])
    succ_gluc_at = e.get_condition(["At"], "Succinate+Glucose", conc * 2, "OD")
    succ_gluc_at.filter_time(xs[-1])
    gluc_at = e.get_condition(["At"], "Glucose", conc, "OD")
    gluc_at.filter_time(xs[-1])
    fig = succ_at.plot_condition()
    at = Species("At", params.loc["At"])
    oa = Species("Oa", params.loc["Oa"])
    oa.N0 = 0
    m = Model(at, oa, e, xs, conc, 0)
    m.integrate_mac_arthur()
    fig.add_trace(
        go.Scatter(
            x=m.t,
            y=m.at.y,
            mode="lines",
            line=dict(width=2, color="blue"),
            name="At Simulation",
        )
    )
    fig = style_plot(fig)
    fig.update_layout(
        title="Succinate At Simulations vs Data",
        xaxis_title="Time (hours)",
        yaxis_title="OD600",
    )
    # fig.write_image("plots/simulations/at_succinate.svg")

    fig = gluc_at.plot_condition()
    at = Species("At", params.loc["At"])
    oa = Species("Oa", params.loc["Oa"])
    oa.N0 = 0
    at.a = 0
    m = Model(at, oa, e, xs, conc, 0)
    m.integrate_mac_arthur()
    fig.add_trace(
        go.Scatter(
            x=m.t,
            y=m.at.y,
            mode="lines",
            line=dict(width=2, color="blue"),
            name="At Simulation",
        )
    )
    fig = style_plot(fig)
    fig.update_layout(
        title="Glucose At Simulations vs Data",
        xaxis_title="Time (hours)",
        yaxis_title="OD600",
    )
    # fig.write_image("plots/simulations/at_glucose.svg")

    fig = succ_gluc_at.plot_condition()
    at = Species("At", params.loc["At"])
    oa = Species("Oa", params.loc["Oa"])
    oa.N0 = 0
    at.a = 1
    m = Model(at, oa, e, xs, conc, 0)
    m.integrate_mac_arthur()
    fig.add_trace(
        go.Scatter(
            x=m.t,
            y=m.at.y,
            mode="lines",
            line=dict(width=2, color="blue"),
            name="At Simulation",
            customdata=np.column_stack([m.succinate, m.glucose]),
            hovertemplate=(
                "t = %{x:.2f} h<br>"
                "At = %{y:.3f} OD<br>"
                "Succinate = %{customdata[0]:.3f} mM<br>"
                "Glucose = %{customdata[1]:.3f} mM<br>"
                "<extra></extra>"
            ),
        )
    )
    fig = style_plot(fig)
    fig.update_layout(
        title="Succinate+Glucose At Simulations vs Data",
        xaxis_title="Time (hours)",
        yaxis_title="OD600",
    )
    fig.show()
    fig.write_image("plots/simulations/at_succinate_glucose.svg")


def fit_oa_glucose_succinate():
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)
    e = Experiment(d="data/251018_succinate_glucose_plate_reader/metaod/")
    xs = np.linspace(0, 28, 100)
    succ_oa = e.get_condition(["Oa"], "Succinate", conc, "OD")
    succ_oa.filter_time(xs[-1])
    succ_gluc_oa = e.get_condition(["Oa"], "Succinate+Glucose", conc * 2, "OD")
    succ_gluc_oa.filter_time(xs[-1])
    gluc_oa = e.get_condition(["Oa"], "Glucose", conc, "OD")
    gluc_oa.filter_time(xs[-1])
    fig = succ_oa.plot_condition()
    at = Species("At", params.loc["At"])
    oa = Species("Oa", params.loc["Oa"])
    at.N0 = 0
    m = Model(at, oa, e, xs, conc, 0)
    m.integrate_mac_arthur()
    fig.add_trace(
        go.Scatter(
            x=m.t,
            y=m.oa.y,
            mode="lines",
            line=dict(width=2, color="blue"),
            name="Oa Simulation",
        )
    )
    fig = style_plot(fig)
    fig.update_layout(
        title="Succinate Oa Simulations vs Data",
        xaxis_title="Time (hours)",
        yaxis_title="OD600",
    )
    # fig.write_image("plots/simulations/oa_succinate.svg")

    fig = gluc_oa.plot_condition()
    at = Species("At", params.loc["At"])
    oa = Species("Oa", params.loc["Oa"])
    at.N0 = 0
    oa.a = 0
    m = Model(at, oa, e, xs, conc, 0)
    m.integrate_mac_arthur()
    fig.add_trace(
        go.Scatter(
            x=m.t,
            y=m.oa.y,
            mode="lines",
            line=dict(width=2, color="blue"),
            name="Oa Simulation",
        )
    )
    fig = style_plot(fig)
    fig.update_layout(
        title="Glucose Oa Simulations vs Data",
        xaxis_title="Time (hours)",
        yaxis_title="OD600",
    )
    # fig.write_image("plots/simulations/oa_glucose.svg")

    fig = succ_gluc_oa.plot_condition()
    at = Species("At", params.loc["At"])
    oa = Species("Oa", params.loc["Oa"])
    at.N0 = 0
    oa.a = 0.75
    m = Model(at, oa, e, xs, conc, 0)
    m.integrate_mac_arthur()
    fig.add_trace(
        go.Scatter(
            x=m.t,
            y=m.oa.y,
            mode="lines",
            line=dict(width=2, color="blue"),
            name="Oa Simulation",
        )
    )
    fig = style_plot(fig)
    fig.update_layout(
        title="Succinate+Glucose Oa Simulations vs Data",
        xaxis_title="Time (hours)",
        yaxis_title="OD600",
    )
    fig.write_image("plots/simulations/oa_succinate_glucose.svg")


fit_oa_glucose_succinate()
