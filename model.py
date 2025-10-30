import pandas as pd
from scipy.integrate import odeint
import plotly.graph_objects as go
from style import style_plot


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

    def simulate(self, y, t):
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
        self.at.y, self.oa.y, self.glucose, self.succinate = odeint(
            self.simulate, y0, self.t
        ).T

    def plot_at_oa(self):
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
