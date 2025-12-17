from scipy.integrate import odeint
from experiment import Species, Experiment
import pandas as pd
from os import path
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from style import style_plot, colors


class Reactor:
    def __init__(self, name, at, oa):
        self.name = name
        self.at = at
        self.oa = oa
        self.glucose = None
        self.succinate = None


class Model:
    def __init__(self, at, oa, e, t, C_mono, D):
        self.c1 = Reactor("C1", at, oa)
        self.c2 = Reactor("C2", at, oa)
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

    def simulate(self, y, t):
        at_c1, oa_c1, g_c1, s_c1, at_c2, oa_c2, g_c2, s_c2 = y
        J_S_at_c1 = (
            self.c1.at.a * self.c1.at.v_succ_lag * s_c1 / (self.c1.at.K_succ + s_c1)
        )
        J_G_at_c1 = (
            (1.0 - self.c1.at.a)
            * self.c1.at.v_gluc_lag
            * g_c1
            / (self.c1.at.K_gluc + g_c1)
        )
        J_S_oa_c1 = (
            self.c1.oa.a * self.c1.oa.v_succ_lag * s_c1 / (self.c1.oa.K_succ + s_c1)
        )
        J_G_oa_c1 = (
            (1.0 - self.c1.oa.a)
            * self.c1.oa.v_gluc_lag
            * g_c1
            / (self.c1.oa.K_gluc + g_c1)
        )

        J_at_c1 = J_S_at_c1 + J_G_at_c1
        J_oa_c1 = J_S_oa_c1 + J_G_oa_c1

        J_S_at_c2 = (
            self.c2.at.a * self.c2.at.v_succ_lag * s_c2 / (self.c2.at.K_succ + s_c2)
        )
        J_G_at_c2 = (
            (1.0 - self.c2.at.a)
            * self.c2.at.v_gluc_lag
            * g_c2
            / (self.c2.at.K_gluc + g_c2)
        )
        J_S_oa_c2 = (
            self.c2.oa.a * self.c2.oa.v_succ_lag * s_c2 / (self.c2.oa.K_succ + s_c2)
        )
        J_G_oa_c2 = (
            (1.0 - self.c2.oa.a)
            * self.c2.oa.v_gluc_lag
            * g_c2
            / (self.c2.oa.K_gluc + g_c2)
        )

        J_at_c2 = J_S_at_c2 + J_G_at_c2
        J_oa_c2 = J_S_oa_c2 + J_G_oa_c2

        datdt_c1 = at_c1 * (J_at_c1 - self.D)
        doadt_c1 = oa_c1 * (J_oa_c1 - self.D)
        datdt_c2 = at_c2 * (J_at_c2 - self.D) + self.D * at_c1
        doadt_c2 = oa_c2 * (J_oa_c2 - self.D) + self.D * oa_c1

        dgdt_c1 = self.D * (self.M_glucose - g_c1) - (
            J_G_at_c1 * at_c1 / self.c1.at.q_gluc
            + J_G_oa_c1 * oa_c1 / self.c1.oa.q_gluc
        )
        dsdt_c1 = self.D * (self.M_succinate - s_c1) - (
            J_S_at_c1 * at_c1 / self.c1.at.q_succ
            + J_S_oa_c1 * oa_c1 / self.c1.oa.q_succ
        )
        dgdt_c2 = self.D * (g_c1 - g_c2) - (
            J_G_at_c2 * at_c2 / self.c2.at.q_gluc
            + J_G_oa_c2 * oa_c2 / self.c2.oa.q_gluc
        )
        dsdt_c2 = self.D * (s_c1 - s_c2) - (
            J_S_at_c2 * at_c2 / self.c2.at.q_succ
            + J_S_oa_c2 * oa_c2 / self.c2.oa.q_succ
        )
        return [
            datdt_c1,
            doadt_c1,
            dgdt_c1,
            dsdt_c1,
            datdt_c2,
            doadt_c2,
            dgdt_c2,
            dsdt_c2,
        ]

    def integrate_model(self):
        y0 = [
            self.c1.at.N0,
            self.c1.oa.N0,
            self.M_glucose,
            self.M_succinate,
            0,
            0,
            self.M_glucose,
            self.M_succinate,
        ]
        sol = odeint(self.simulate, y0, self.t)

        self.c1.at_y, self.c1.oa_y, self.c1.glucose, self.c1.succinate = (
            sol[:, 0],
            sol[:, 1],
            sol[:, 2],
            sol[:, 3],
        )
        self.c2.at_y, self.c2.oa_y, self.c2.glucose, self.c2.succinate = (
            sol[:, 4],
            sol[:, 5],
            sol[:, 6],
            sol[:, 7],
        )


xs = np.linspace(0, 100, 1000)
p_f = path.join("parameters", f"parameters_{15}_mM_C.csv")
params = pd.read_csv(p_f, index_col=0)
at = Species("At", params.loc["At"])
at.N0 = 0.221667
at.a = 0.9
oa = Species("Oa", params.loc["Oa"])
oa.N0 = 0
m = Model(at, oa, None, xs, 15, 0.3)
m.integrate_model()

fig = make_subplots(
    rows=2, cols=1, shared_xaxes=True, subplot_titles=("C1", "C2"), shared_yaxes=True
)

for i, c in enumerate([m.c1, m.c2]):
    fig.add_trace(
        go.Scatter(
            x=m.t,
            y=c.at_y,
            mode="lines",
            line=dict(width=2, color=colors["Succinate+Glucose"]),
            name=f"Model At {c.name}",
            showlegend=False,
        ),
        row=i + 1,
        col=1,
    )
fig.update_layout(
    xaxis=dict(title="Time [h]", ticks="inside"),
    yaxis=dict(title="OD600", ticks="inside"),
)
fig = style_plot(fig, line_thickness=2)
fig.write_image("tmp.svg")
