from sympy import symbols, Eq, solve
from sympy import init_printing
import pandas as pd
from os import path
from experiment import Species

init_printing(use_unicode=True)
p_f = path.join("parameters", f"parameters_{15}_mM_C.csv")
params = pd.read_csv(p_f, index_col=0)
at = Species("At", params.loc["At"])
oa = Species("Oa", params.loc["Oa"])
D = 0.3

(
    at_c1,
    oa_c1,
    s_c1,
    g_c1,
    at_c2,
    oa_c2,
    s_c2,
    g_c2,
    J_at_c1,
    J_oa_c1,
    J_at_c2,
    J_oa_c2,
) = symbols(
    "at_c1 oa_c1 s_c1 g_c1 at_c2 oa_c2 s_c2 g_c2 J_at_c1 J_oa_c1 J_at_c2 J_oa_c2"
)

eq_J_at_c1 = Eq(
    J_at_c1,
    at.a * at.v_gluc * g_c1 / (at.K_gluc + g_c1)
    + (1 - at.a) * at.v_succ * s_c1 / (at.K_succ + s_c1),
)
eq_J_oa_c1 = Eq(
    J_oa_c1,
    oa.a * oa.v_gluc * g_c1 / (oa.K_gluc + g_c1)
    + (1 - oa.a) * oa.v_succ * s_c1 / (oa.K_succ + s_c1),
)
eq_J_at_c2 = Eq(
    J_at_c2,
    at.a * at.v_gluc * g_c2 / (at.K_gluc + g_c2)
    + (1 - at.a) * at.v_succ * s_c2 / (at.K_succ + s_c2),
)
eq_J_oa_c2 = Eq(
    J_oa_c2,
    oa.a * oa.v_gluc * g_c2 / (oa.K_gluc + g_c2)
    + (1 - oa.a) * oa.v_succ * s_c2 / (oa.K_succ + s_c2),
)
