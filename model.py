import pandas as pd
from scipy.integrate import odeint
import plotly.graph_objects as go
from style import style_plot
import numpy as np
from scipy.integrate import solve_ivp


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

        self.M_succinate = self.C_to_mM_succinate[self.C_mono]
        self.M_glucose = self.C_to_mM_glucose[self.C_mono]
        self.glucose = None
        self.succinate = None

    def simulate(self, t, y):
        at, oa, g, s = y
        J_S_at = self.at.a * self.at.v_succ_lag * s / (self.at.K_succ + s)
        J_G_at = (1.0 - self.at.a) * self.at.v_gluc_lag * g / (self.at.K_gluc + g)
        J_S_oa = self.oa.a * self.oa.v_succ_lag * s / (self.oa.K_succ + s)
        J_G_oa = (1.0 - self.oa.a) * self.oa.v_gluc_lag * g / (self.oa.K_gluc + g)

        J_at = J_S_at + J_G_at
        J_oa = J_S_oa + J_G_oa

        # D may be zero in batch; keep form general
        datdt = at * (J_at - self.D)
        doadt = oa * (J_oa - self.D)

        dgdt = self.D * (self.M_glucose - g) - (
            J_G_at * at / self.at.q_gluc + J_G_oa * oa / self.oa.q_gluc
        )
        dsdt = self.D * (self.M_succinate - s) - (
            J_S_at * at / self.at.q_succ + J_S_oa * oa / self.oa.q_succ
        )

        return [datdt, doadt, dgdt, dsdt]

    def integrate_model(self):
        y0 = [self.at.N0, self.oa.N0, self.M_glucose, self.M_succinate]
        sol = solve_ivp(
            self.simulate,
            (self.t[0], self.t[-1]),
            y0,
            t_eval=self.t,
            method="BDF",
            rtol=1e-8,
            atol=1e-12,
        )

        Y = sol.y.T
        Y = sol.y.T
        Y[np.abs(Y) < 1e-12] = 0.0
        self.at.y, self.oa.y, self.glucose, self.succinate = (
            Y[:, 0],
            Y[:, 1],
            Y[:, 2],
            Y[:, 3],
        )

    def plot_at_oa(self):
        fig = go.Figure()
        species = []
        if (self.M_glucose != 0) and (self.M_succinate != 0):
            customdata = np.column_stack([self.succinate, self.glucose])
            hovertemplate = (
                "Time: %{x:.2f} hours<br>"
                "Succinate: %{customdata[0]:.2f} mM<br>"
                "Glucose: %{customdata[1]:.2f} mM<br>"
                "<extra></extra>"
            )
        elif (self.M_glucose != 0) and (self.M_succinate == 0):
            customdata = self.glucose
            hovertemplate = (
                "Time: %{x:.2f} hours<br>"
                "Glucose: %{customdata:.2f} mM<br>"
                "<extra></extra>"
            )
        elif (self.M_glucose == 0) and (self.M_succinate != 0):
            customdata = self.succinate
            hovertemplate = (
                "Time: %{x:.2f} hours<br>"
                "Succinate: %{customdata:.2f} mM<br>"
                "<extra></extra>"
            )
        if self.at.y[0] != 0:
            fig.add_trace(
                go.Scatter(
                    x=self.t,
                    y=self.at.y,
                    mode="lines",
                    line=dict(width=2, color="green"),
                    name="Model At",
                    customdata=customdata,
                    hovertemplate=hovertemplate,
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
                    customdata=customdata,
                    hovertemplate=hovertemplate,
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
                    customdata=customdata,
                    hovertemplate=hovertemplate,
                )
            )
        title = f"{' + '.join(species)} - {self.C_to_mM_glucose[self.C_mono]} mM Glucose, {self.C_to_mM_succinate[self.C_mono]} mM Succinate"
        fig.update_layout(
            title=title,
            xaxis_title="Time (hours)",
            yaxis_title="OD600",
            width=600,
            height=400,
        )
        fig = style_plot(fig, font_size=14, line_thickness=2)
        return fig

    def plot_substrates(self):
        fig = go.Figure()
        if self.glucose[0] != 0:
            fig.add_trace(
                go.Scatter(
                    x=self.t,
                    y=self.glucose,
                    mode="lines",
                    line=dict(width=2, color="blue"),
                    name="Model Glucose",
                )
            )
        if self.succinate[0] != 0:
            fig.add_trace(
                go.Scatter(
                    x=self.t,
                    y=self.succinate,
                    mode="lines",
                    line=dict(width=2, color="orange"),
                    name="Model Succinate",
                )
            )
        title = f"Substrates - {self.C_to_mM_glucose[self.C_mono]} mM Glucose, {self.C_to_mM_succinate[self.C_mono]} mM Succinate"
        fig.update_layout(
            title=title,
            xaxis_title="Time (hours)",
            yaxis_title="Concentration (mM)",
            width=600,
            height=400,
        )
        fig = style_plot(fig, font_size=14, line_thickness=2)
        return fig


"""
from os import path
from experiment import Species, Experiment

xs = np.linspace(0, 100, 1000)
p_f = path.join("parameters", f"parameters_{15}_mM_C.csv")
params = pd.read_csv(p_f, index_col=0)
at = Species("At", params.loc["At"])
at.N0 = 0.221667
at.a = 0.94736842
oa = Species("Oa", params.loc["Oa"])
oa.N0 = 0.0
# oa.a = 0.0

m = Model(at, oa, None, xs, 15, 0.3)
m.integrate_model()
fig = m.plot_at_oa()
fig.write_image("tmp.svg")
"""
