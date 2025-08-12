from sympy import symbols, Eq, solve
from sympy import init_printing
import pandas as pd

init_printing(use_unicode=True)
df = dict(pd.read_csv("parameters.csv"))
p = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
R1C1, R1C2, N1C1, N1C2, N1JC1, N1JC2, dN1C1, dN1C2, dR1C1, dR1C2 = symbols(
    "R1C1 R1C2 N1C1 N1C2 N1JC1 N1JC2 dN1C1 dN1C2 dR1C1 dR1C2"
)

eq_N1JC1 = Eq(N1JC1, p["v1_1"] * R1C1 / (R1C1 + p["K1_1"]))
eq_N1JC2 = Eq(N1JC2, p["v1_1"] * R1C2 / (R1C2 + p["K1_1"]))
eq_N1C1 = Eq(dN1C1, N1C1 * (eq_N1JC1.rhs - p["D"]))
eq_N1C2 = Eq(dN1C2, N1C2 * (eq_N1JC2.rhs - p["D"]) + p["D"] * N1C1)
eq_R1C1 = Eq(dR1C1, p["D"] * (p["M1"] - R1C1) - eq_N1JC1.rhs * N1C1 / p["q1_1"])
eq_R1C2 = Eq(dR1C2, p["D"] * (R1C1 - R1C2) - eq_N1JC2.rhs * N1C2 / p["q1_1"])

R1fC1 = solve(Eq(p["D"], eq_N1JC1.rhs), R1C1)[0]
N1fC1 = solve(Eq(0, eq_R1C1.rhs), N1C1)[0].subs({R1C1: R1fC1})
eq_D_eff = Eq(eq_N1JC2.rhs, p["D"] * (1 - N1C1 / N1C2))


def solve_D_eff():
    sol1, sol2 = solve(
        Eq(eq_D_eff.lhs.subs({R1C2: solve(Eq(0, eq_R1C2.rhs), R1C2)[0]}), eq_D_eff.rhs),
        N1C2,
    )
    sol1 = sol1.subs({R1C1: R1fC1, N1C1: N1fC1})
    sol2 = sol2.subs({R1C1: R1fC1, N1C1: N1fC1})
    sol = sol1 if sol1 > 0 else sol2
    D_eff = p["D"] * (1 - N1fC1 / sol)
    R1fC2 = solve(Eq(D_eff, eq_N1JC2.rhs), R1C2)[0]
    N1fC2 = solve(Eq(0, eq_R1C2.rhs), N1C2)[0].subs({R1C2: R1fC2, R1C1: R1fC1})
    print("N1C1:", N1fC1, "N1C2:", N1fC2)
    print("R1C1:", R1fC1, "R1C2:", R1fC2)
