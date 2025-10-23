from experiment import Condition, Experiment, Species
import pandas as pd
from scipy.integrate import odeint
import plotly.graph_objects as go
from style import style_plot
import numpy as np
from os import path


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

    def lag_time(self, t, t_lag, lag_width):
        if t_lag <= 0:
            return 1.0
        k = 4.394 / max(lag_width, 1e-9)
        return 1.0 / (1.0 + np.exp(-k * (t - t_lag)))

    def mac_arthur_lag(self, y, t):
        at, oa, g, s = y

        # ON/OFF gates
        phi_S_at = (
            1
            if self.at.lag_succ == 0
            else self.lag_time(t, self.at.lag_succ, self.at.lag_succ_width)
        )
        phi_G_at = (
            1
            if self.at.lag_gluc == 0
            else self.lag_time(t, self.at.lag_gluc, self.at.lag_gluc_width)
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

        # Uptake fluxes (succinate uses a; glucose uses 1-a)
        J_S_at = phi_S_at * self.at.a * self.at.v_succ * s / (self.at.K_succ + s)
        J_G_at = (
            phi_G_at * (1.0 - self.at.a) * self.at.v_gluc * g / (self.at.K_gluc + g)
        )
        J_S_oa = phi_S_oa * self.oa.a * self.oa.v_succ * s / (self.oa.K_succ + s)
        J_G_oa = (
            phi_G_oa * (1.0 - self.oa.a) * self.oa.v_gluc * g / (self.oa.K_gluc + g)
        )

        J_at = J_S_at + J_G_at
        J_oa = J_S_oa + J_G_oa
        print(J_at)

        datdt = at * (J_at - self.D)
        doadt = oa * (J_oa - self.D)

        dgdt = self.D * (self.M_glucose - g) - (
            J_G_at * at / self.at.q_gluc + J_G_oa * oa / self.oa.q_gluc
        )
        dsdt = self.D * (self.M_succinate - s) - (
            J_S_at * at / self.at.q_succ + J_S_oa * oa / self.oa.q_succ
        )

        return [datdt, doadt, dgdt, dsdt]

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
                    line=dict(width=2, color="purple"),
                    name="Model",
                )
            )
            species.append("At")
        if self.oa.y[0] != 0:
            fig.add_trace(
                go.Scatter(
                    x=self.t,
                    y=self.oa.y,
                    mode="lines",
                    line=dict(width=2, color="purple"),
                    name="Model",
                )
            )
            species.append("Oa")
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


def plot_fit():
    e = Experiment(d="data/251018_succinate_glucose_plate_reader/metaod/")
    concentrations = [2.5, 5, 7.5, 15]
    carbon_sources = ["Succinate", "Glucose"]

    for conc in concentrations:
        # load parameter file once per concentration
        p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
        params = pd.read_csv(p_f, index_col=0)

        for cs in carbon_sources:
            a_pref = 1.0 if cs == "Succinate" else 0.0

            for focal in ("At",):
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
                tgrid = cond.xs[0]
                print(conc)
                M = Model(at, oa, e, tgrid, conc, 0)
                M.integrate_mac_arthur()
                M.plot_at_oa(fig)

                # pick the correct kinetic params for the title
                sp = at if focal == "At" else oa
                if cs == "Succinate":
                    mu_val, K_val = sp.v_succ, sp.K_succ
                else:
                    mu_val, K_val = sp.v_gluc, sp.K_gluc

                title = (
                    f"{focal} on {cs} {conc:g} mM — μ: {mu_val:g} 1/h   K: {K_val:g} mM"
                )
                fig.update_layout(title=title)
                fig = style_plot(fig, font_size=14, line_thickness=2)

                safe_conc = f"{conc:g}".replace(".", "_")
                outfile = path.join(
                    "plots",
                    "fitting",
                    f"{focal.lower()}_{safe_conc}_mM_{cs.lower()}.svg",
                )
                fig.write_image(outfile)


plot_fit()
