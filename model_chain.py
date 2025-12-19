from scipy.integrate import odeint
from experiment import Species, Experiment
import pandas as pd
from os import path
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from style import style_plot, colors
from scipy.integrate import solve_ivp


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

    def simulate(self, t, y):
        at_c1, oa_c1, g_c1, s_c1, at_c2, oa_c2, g_c2, s_c2 = y

        epsN = 1e-12
        g_c1 = max(g_c1, 0.0)
        s_c1 = max(s_c1, 0.0)
        g_c2 = max(g_c2, 0.0)
        s_c2 = max(s_c2, 0.0)
        at_c1 = 0.0 if at_c1 < epsN else at_c1
        oa_c1 = 0.0 if oa_c1 < epsN else oa_c1
        at_c2 = 0.0 if at_c2 < epsN else at_c2
        oa_c2 = 0.0 if oa_c2 < epsN else oa_c2

        def monod(s, K):
            if s <= 0.0:
                return 0.0
            return s / (K + s)

        fS_at_c1 = monod(s_c1, self.c1.at.K_succ)
        fG_at_c1 = monod(g_c1, self.c1.at.K_gluc)
        fS_oa_c1 = monod(s_c1, self.c1.oa.K_succ)
        fG_oa_c1 = monod(g_c1, self.c1.oa.K_gluc)

        J_S_at_c1 = self.c1.at.a * self.c1.at.v_succ_lag * fS_at_c1
        J_G_at_c1 = (1.0 - self.c1.at.a) * self.c1.at.v_gluc_lag * fG_at_c1
        J_S_oa_c1 = self.c1.oa.a * self.c1.oa.v_succ_lag * fS_oa_c1
        J_G_oa_c1 = (1.0 - self.c1.oa.a) * self.c1.oa.v_gluc_lag * fG_oa_c1

        J_at_c1 = J_S_at_c1 + J_G_at_c1
        J_oa_c1 = J_S_oa_c1 + J_G_oa_c1

        fS_at_c2 = monod(s_c2, self.c2.at.K_succ)
        fG_at_c2 = monod(g_c2, self.c2.at.K_gluc)
        fS_oa_c2 = monod(s_c2, self.c2.oa.K_succ)
        fG_oa_c2 = monod(g_c2, self.c2.oa.K_gluc)

        J_S_at_c2 = self.c2.at.a * self.c2.at.v_succ_lag * fS_at_c2
        J_G_at_c2 = (1.0 - self.c2.at.a) * self.c2.at.v_gluc_lag * fG_at_c2
        J_S_oa_c2 = self.c2.oa.a * self.c2.oa.v_succ_lag * fS_oa_c2
        J_G_oa_c2 = (1.0 - self.c2.oa.a) * self.c2.oa.v_gluc_lag * fG_oa_c2

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
            0.0,
            0.0,
            self.M_glucose,
            self.M_succinate,
        ]

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
        Y[~np.isfinite(Y)] = np.nan
        Y[Y < 0] = 0.0
        Y[np.abs(Y) < 1e-12] = 0.0

        self.c1.at_y, self.c1.oa_y, self.c1.glucose, self.c1.succinate = (
            Y[:, 0],
            Y[:, 1],
            Y[:, 2],
            Y[:, 3],
        )
        self.c2.at_y, self.c2.oa_y, self.c2.glucose, self.c2.succinate = (
            Y[:, 4],
            Y[:, 5],
            Y[:, 6],
            Y[:, 7],
        )

    def plot_at_oa(self):
        fig = make_subplots(
            rows=2,
            cols=1,
            subplot_titles=("C1", "C2"),
            shared_xaxes=True,
        )
        species = []
        for i, c in enumerate([self.c1, self.c2]):
            print(i)
            if (self.M_glucose != 0) and (self.M_succinate != 0):
                customdata = np.column_stack([c.succinate, c.glucose])
                hovertemplate = (
                    "Time: %{x:.2f} hours<br>"
                    "Succinate: %{customdata[0]:.2f} mM<br>"
                    "Glucose: %{customdata[1]:.2f} mM<br>"
                    "<extra></extra>"
                )
            elif (self.M_glucose != 0) and (self.M_succinate == 0):
                customdata = c.glucose
                hovertemplate = (
                    "Time: %{x:.2f} hours<br>"
                    "Glucose: %{customdata:.2f} mM<br>"
                    "<extra></extra>"
                )
            elif (self.M_glucose == 0) and (self.M_succinate != 0):
                customdata = c.succinate
                hovertemplate = (
                    "Time: %{x:.2f} hours<br>"
                    "Succinate: %{customdata:.2f} mM<br>"
                    "<extra></extra>"
                )
            if c.at_y[-1] != 0:
                fig.add_trace(
                    go.Scatter(
                        x=self.t,
                        y=c.at_y,
                        mode="lines",
                        line=dict(width=2, color=colors["at"]),
                        name="Model At",
                        customdata=customdata,
                        hovertemplate=hovertemplate,
                        showlegend=False,
                    ),
                    row=i + 1,
                    col=1,
                )
                species.append("At")
            if c.oa_y[-1] != 0:
                fig.add_trace(
                    go.Scatter(
                        x=self.t,
                        y=c.oa_y,
                        mode="lines",
                        line=dict(width=2, color=colors["oa"]),
                        name="Model Oa",
                        customdata=customdata,
                        hovertemplate=hovertemplate,
                        showlegend=False,
                    ),
                    row=i + 1,
                    col=1,
                )
                species.append("Oa")
            if c.oa_y[-1] != 0 and c.at_y[-1] != 0:
                fig.add_trace(
                    go.Scatter(
                        x=self.t,
                        y=c.at_y + c.oa_y,
                        mode="lines",
                        line=dict(width=2, color="purple"),
                        name="Model Total",
                        customdata=customdata,
                        hovertemplate=hovertemplate,
                        showlegend=False,
                    ),
                    row=i + 1,
                    col=1,
                )
        legend = set(species)
        if len(legend) == 2:
            legend.add("Total")
        for l in legend:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="lines",
                    name=l,
                    marker=dict(color=colors[l]),
                )
            )

        fig.update_layout(
            xaxis_title="Time (hours)",
            yaxis_title="OD600",
            width=600,
            height=400,
        )
        fig = style_plot(fig, font_size=14, line_thickness=2)
        return fig


"""xs = np.linspace(0, 100, 1000)
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
