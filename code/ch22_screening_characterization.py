"""Original screening and diffraction audits; synthetic phases, not an \
assessed alloy."""

from pathlib import Path
from itertools import product
from math import comb, exp, log, sin, cos, asin, sqrt, pi
import numpy as np
import pycalphad
from pycalphad import Database, equilibrium, variables as v

assert pycalphad.__version__ == "0.11.1"
R = float(v.R)
# Integer atomic-fraction grid: all five components present, 5 through 50 at.%.
grid = [q for q in product(range(1, 11), repeat=5) if sum(q) == 20]
with_zero = sum(sum(q) == 20 for q in product(range(11), repeat=5))
assert len(grid) == comb(19, 4) - 5 * comb(9, 4) == 3246
assert with_zero == comb(24, 4) - 5 * comb(13, 4) == 7051
print(
    (
        f"Five-component 5%-grid, positive / zero-allowed: "
        f"{len(grid)}"
        f", "
        f"{with_zero}"
    )
)
print(
    (
        f"Unrestricted 26-element equimolar subsets of sizes3--6: "
        f"{sum((comb(26, k) for k in range(3, 7)))}"
    )
)
tm_c = 1250.0
print(
    (
        f"0.8 times 1250 Celsius: operational="
        f"{0.8 * tm_c:.3f}"
        f" C; absolute-ratio="
        f"{0.8 * (tm_c + 273.15) - 273.15:.3f}"
        f" C"
    )
)
print(
    (
        f"Endpoint-only fraction test: endpoints=.100,.100; interior="
        f"{0.1 + 0.6 * sin(pi / 2) ** 2:.3f}"
    )
)

# An infinitesimal product can partition while the remaining matrix conserves
# composition.
A, z, xp, advantage = 10000.0, 0.3, 0.7, 300.0
same_composition = A * (z - xp) ** 2 - advantage
print(
    (
        f"Parent tangent-plane advantage="
        f"{advantage:.6f}"
        f"; same-composition product penalty="
        f"{same_composition:.6f}"
        f" J/mol"
    )
)
for amount in (0.1, 0.2):
    xm = (z - amount * xp) / (1 - amount)
    energy = (1 - amount) * A * (xm - z) ** 2 - amount * advantage
    assert abs((1 - amount) * xm + amount * xp - z) < 1e-14
    print(
        (
            f"Partitioned product amount="
            f"{amount:.1f}"
            f": matrix x="
            f"{xm:.9f}"
            f", total change="
            f"{energy:.9f}"
            f" J/mol"
        )
    )

# Execute a small batch with two ideal substitutional solution branches.
db = Database(str(Path("data/ch22_screening.tdb")))
D = 6000.0
temperatures = (900.0, 1000.0)
compositions = [j / 20 for j in range(7, 14)]
rows = {x: [] for x in compositions}
for T in temperatures:
    delta = 8 * (1100 - T)
    u, w = exp((D + delta) / (R * T)), exp((delta - D) / (R * T))
    b = (u - 1) / (u - w)
    a = w * b
    mu_a = R * T * log(1 - a)
    slope = D + R * T * log(a / (1 - a))
    print(
        (
            f"Synthetic common-line contacts T="
            f"{T:.0f}"
            f" K: MATRIX="
            f"{a:.9f}"
            f", PRECIP="
            f"{b:.9f}"
        )
    )
    for x in compositions:
        result = equilibrium(
            db,
            ["CU", "NI"],
            ["MATRIX", "PRECIP"],
            {v.T: T, v.P: 100000, v.N: 1, v.X("NI"): x},
            calc_opts={"pdens": 200},
        )
        records = [
            (str(p), float(c), float(n))
            for p, c, n in zip(
                result.Phase.values.reshape(-1),
                result.X.sel(component="NI").values.reshape(-1),
                result.NP.values.reshape(-1),
            )
            if p and np.isfinite(n) and n > 1e-9
        ]
        assert abs(sum(n for _, _, n in records) - 1) < 1e-9
        assert abs(sum(c * n for _, c, n in records) - x) < 1e-9
        fraction = sum(n for p, _, n in records if p == "PRECIP")
        analytic = max(0.0, min(1.0, (x - a) / (b - a)))
        assert abs(fraction - analytic) < 1e-7
        if a < x < b:
            assert (
                abs(float(result.GM.values.squeeze()) - (mu_a + slope * x))
                < 1e-6
            )
            assert all(
                abs(c - ({"MATRIX": a, "PRECIP": b}[p])) < 1e-7
                for p, c, _ in records
            )
        rows[x].append(fraction)
nominal = [x for x, fs in rows.items() if all(0.05 < f < 0.5 for f in fs)]
# A stipulated +/-0.04 bound on each fraction; NOT an inferred uncertainty.
robust = [
    x
    for x, fs in rows.items()
    if all(0.05 < f - 0.04 and f + 0.04 < 0.5 for f in fs)
]
for x, fs in rows.items():
    print(
        (
            f"Synthetic xNi="
            f"{x:.2f}"
            f": fPRECIP900/1000="
            f"{fs[0]:.6f}"
            f","
            f"{fs[1]:.6f}"
            f"; nominal="
            f"{x in nominal}"
            f"; robust="
            f"{x in robust}"
        )
    )
assert nominal == [0.5, 0.55] and robust == [0.5]
print(
    (
        "Executed equilibrium batch:14 cases, every mass/composition balance "
        "and analytic fraction PASS"
    )
)


# Classifier behavior for an absent calculation, not a claim of an actual
# solver failure.
def classify(values):
    if values is None:
        return "unknown"
    return "pass" if all(0.05 < f < 0.5 for f in values) else "reject"


assert classify(None) == "unknown"
print("Deliberately absent result classified as unknown, not reject or pass")

# Synthetic angular diffraction data; known cubic indexing and wavelength.
wavelength, a0, true_a, zero = (
    0.15406,
    0.36,
    0.3618,
    0.04,
)  # nm, nm, nm, degrees
squares = (3, 4, 8, 11, 12)


def angles(a):
    return [
        2 * asin(wavelength * sqrt(q) / (2 * a)) * 180 / pi for q in squares
    ]


old = angles(a0)
observed = [t + zero for t in angles(true_a)]
for q, old_t, new_t in zip(squares[:3], old[:3], angles(true_a)[:3]):
    linear = (
        -2 * sin(old_t * pi / 360) / cos(old_t * pi / 360) * (0.005) * 180 / pi
    )
    print(
        (
            f"Cubic h2+k2+l2="
            f"{q}"
            f": exact 2theta shift="
            f"{new_t - old_t:.6f}"
            f" deg; first-order="
            f"{linear:.6f}"
            f" deg"
        )
    )


# For each candidate a, the best common zero offset is the mean angular
# residual.
def fit_at(a):
    residual = [o - t for o, t in zip(observed, angles(a))]
    shift = sum(residual) / len(residual)
    return sum((r - shift) ** 2 for r in residual), shift


left, right = 0.355, 0.367
gold = (sqrt(5) - 1) / 2
c = right - gold * (right - left)
d = left + gold * (right - left)
for _ in range(100):
    if fit_at(c)[0] < fit_at(d)[0]:
        right, d = d, c
        c = right - gold * (right - left)
    else:
        left, c = c, d
        d = left + gold * (right - left)
fitted = (left + right) / 2
sse, fitted_zero = fit_at(fitted)
assert abs(fitted - true_a) < 1e-11 and abs(fitted_zero - zero) < 1e-8
assert sse < 1e-20
print(
    (
        f"Synthetic five-peak joint fit: a="
        f"{fitted:.9f}"
        f" nm; zero="
        f"{fitted_zero:.9f}"
        f" deg; SSE<1e-20 deg^2"
    )
)
print(
    (
        "Scope: no measured pattern, real-alloy screening, mechanical-propert"
        "y prediction or thermal-history simulation"
    )
)
