"""Actually executed pycalphad teaching model; NOT assessed Cu-Ni data.\
\nRun from the book root after installing requirements.txt in a proj\
ect environment.\n"""

from pathlib import Path
from math import exp, log, isclose
import numpy as np
import pycalphad
from pycalphad import Database, calculate, equilibrium, variables as v

DB_PATH = Path("data/ch14_synthetic.tdb")
db = Database(str(DB_PATH))
R = float(v.R)  # Match this engine's rounded constant for the analytic audit.
assert pycalphad.__version__ == "0.11.1"
print(f"Engine: pycalphad {pycalphad.__version__}; R={R:.4f} J/(mol K)")

# CIAAW tabulated central values used as fixed illustrative molar masses.
Mfe, Mcr = 55.845, 51.9961
wcr = 0.45
xcr = (wcr / Mcr) / ((1 - wcr) / Mfe + wcr / Mcr)
Mbar = 1 / ((1 - wcr) / Mfe + wcr / Mcr)
print(
    f"Fe-45 mass percent Cr: xCr={xcr:.9f}, mean molar mass={Mbar:.9f} g/mol"
)
fmol = [0.75, 0.25]
Mphase = [60.0, 90.0]
Vphase = [1e-5, 2e-5]
mass_beta = fmol[1] * Mphase[1] / sum(f * m for f, m in zip(fmol, Mphase))
volume_beta = fmol[1] * Vphase[1] / sum(f * vv for f, vv in zip(fmol, Vphase))
print(
    (
        f"Teaching fractions: beta mole="
        f"{fmol[1]:.6f}"
        f", mass="
        f"{mass_beta:.6f}"
        f", volume="
        f"{volume_beta:.6f}"
    )
)


def solve(phases, T=600.0, z=0.4, pdens=1000):
    return equilibrium(
        db,
        ["CU", "NI"],
        phases,
        {v.T: T, v.P: 100000, v.N: 1, v.X("NI"): z},
        calc_opts={"pdens": pdens},
    )


def occupied(eq):
    names = eq.Phase.values.reshape(-1)
    amounts = eq.NP.values.reshape(-1)
    xs = eq.X.sel(component="NI").values.reshape(-1)
    return sorted(
        (str(p), float(x), float(f))
        for p, x, f in zip(names, xs, amounts)
        if p and np.isfinite(f) and f > 1e-9
    )


T, z = 600.0, 0.4
eq = solve(["ALPHA", "BETA"], T, z)
rows = occupied(eq)
a = 1 / (1 + exp(6000 / (R * T)))
b = 1 - a
fb = (z - a) / (b - a)
assert [r[0] for r in rows] == ["ALPHA", "BETA"]
assert abs(rows[0][1] - a) < 1e-7 and abs(rows[1][1] - b) < 1e-7
assert abs(rows[1][2] - fb) < 1e-7
assert abs(sum(f for _, _, f in rows) - 1) < 1e-9
assert abs(sum(x * f for _, x, f in rows) - z) < 1e-9
print(
    (
        f"Restricted 600 K: xAlpha="
        f"{rows[0][1]:.6f}"
        f", xBeta="
        f"{rows[1][1]:.6f}"
        f", fBeta="
        f"{rows[1][2]:.6f}"
        f", GM="
        f"{float(eq.GM.values.squeeze()):.4f}"
    )
)
assert abs(float(eq.GM.values.squeeze()) - R * T * log(1 - a)) < 1e-6
print("Analytic coexistence and elemental/amount balances: PASS")

# calculate() samples phase states; it does not perform global equilibrium.
calc = calculate(
    db, ["CU", "NI"], ["ALPHA"], T=T, P=100000, points=np.array([[0.6, 0.4]])
)
g_alpha = float(calc.GM.values.squeeze())
assert (
    abs(g_alpha - (6000 * z + R * T * ((1 - z) * log(1 - z) + z * log(z))))
    < 1e-6
)
print(
    (
        f"Homogeneous ALPHA at xNi=.4: GM="
        f"{g_alpha:.4f}"
        f"; excess above restricted equilibrium="
        f"{g_alpha - float(eq.GM.values.squeeze()):.4f}"
    )
)

all_eq = solve(["ALPHA", "BETA", "GAMMA"])
all_rows = occupied(all_eq)
assert len(all_rows) == 1 and all_rows[0][0] == "GAMMA"
assert (
    abs(
        float(all_eq.GM.values.squeeze())
        - (-5000 + R * T * ((1 - z) * log(1 - z) + z * log(z)))
    )
    < 1e-6
)
assert float(all_eq.GM.values.squeeze()) <= float(eq.GM.values.squeeze())
print(
    (
        f"GAMMA entered: fraction="
        f"{all_rows[0][2]:.6f}"
        f", GM="
        f"{float(all_eq.GM.values.squeeze()):.4f}"
    )
)

for temperature in [400.0, 600.0, 1000.0]:
    sweep = occupied(solve(["ALPHA", "BETA"], temperature))
    fraction = sum(f for p, _, f in sweep if p == "BETA")
    print(f"Fixed bulk sweep T={temperature:.0f} K: fBeta={fraction:.6f}")

# One model, two composition sets: stable DELTA demixing below its critical T.
delta = solve(["DELTA"])
dr = occupied(delta)
assert len(dr) == 2 and all(p == "DELTA" for p, _, _ in dr)
assert abs(sum(x * f for _, x, f in dr) - z) < 1e-8
# Independent bracketed analytic binodal with the same R.
lo, hi = 0.01, 0.3
for _ in range(80):
    mid = (lo + hi) / 2
    value = R * T * log(mid / (1 - mid)) + 12000 * (1 - 2 * mid)
    if value > 0:
        hi = mid
    else:
        lo = mid
binodal = (lo + hi) / 2
assert abs(dr[0][1] - binodal) < 1e-6
print(
    (
        f"DELTA composition sets: xLow="
        f"{dr[0][1]:.6f}"
        f", xHigh="
        f"{dr[1][1]:.6f}"
        f", fHigh="
        f"{dr[1][2]:.6f}"
    )
)

# Refinement is a numerical check; the analytic lower bound supplies
# certification.
refined = solve(["ALPHA", "BETA"], pdens=2000)
assert (
    abs(float(refined.GM.values.squeeze()) - float(eq.GM.values.squeeze()))
    < 1e-6
)
print("Sampling-density repeat against analytic certificate: PASS")

# Off-section tie line: constant bulk composition is not a per-phase
# constraint.
alpha = np.array([0.7, 0.2, 0.1])
beta = np.array([0.1, 0.5, 0.4])
f = 0.4
bulk = (1 - f) * alpha + f * beta
assert np.allclose(bulk, [0.46, 0.32, 0.22])
assert alpha[2] != bulk[2] and beta[2] != bulk[2]
print(
    (
        "Ternary off-section balance: bulk=(.46,.32,.22), phase C fractions=."
        "10 and .40"
    )
)
