"""Synthetic CALPHAD building blocks; no assessed elemental data."""

from math import log, log1p, isclose, comb

R = 8.31446261815324
# Numerical coefficients assume kelvin, J/mol and the declared log scale.
a, b, c, d2, dm1 = -10000.0, 20.0, -25.0, -0.002, 100000.0


def unary(T):
    g = a + b * T + c * T * log(T / 1000) + d2 * T * T + dm1 / T
    s = -b - c * (log(T / 1000) + 1) - 2 * d2 * T + dm1 / T**2
    h = a - c * T - d2 * T * T + 2 * dm1 / T
    cp = -c - 2 * d2 * T - 2 * dm1 / T**2
    return g, s, h, cp


T = 900.0
g, s, h, cp = unary(T)
assert isclose(h, g + T * s, abs_tol=1e-8)
print(f"Unary at 900 K: g={g:.9f}, s={s:.9f}, h={h:.9f}, cp={cp:.9f}")
# A slope shift changes entropy, not enthalpy or heat capacity.
delta_b = 3.0
assert isclose((g + delta_b * T) + T * (s - delta_b), h, abs_tol=1e-8)
print(
    (
        f"Added 3T: delta g="
        f"{delta_b * T:.6f}"
        f", delta s="
        f"{-delta_b:.6f}"
        f", delta h=0, delta cp=0"
    )
)
# A deliberately bad upper branch with an extra 100 J/mol offset.
print(
    (
        "Piecewise join: added 100 J/mol upper offset causes delta g=100, del"
        "ta h=100; not a valid smooth join"
    )
)

x = 0.2
y = 1 - 2 * x
w = x * (1 - x)
L0 = 12000 - 4 * T + 2 * T * log(T / 1000)
L1 = -3000 + T
L2 = 2000
L0p = -2 + 2 * log(T / 1000)
L1p = 1.0
ge = w * (L0 + L1 * y + L2 * y * y)
se = -w * (L0p + L1p * y)
he = ge + T * se
cpe = -2 * w
print(
    f"RK at 900 K, xB=.2: gE={ge:.9f}, sE={se:.9f}, hE={he:.9f}, cpE={cpe:.9f}"
)


# Relabeling A/B reverses odd coefficients; retaining them changes the model.
def rk(x, coeff):
    return x * (1 - x) * sum(L * (1 - 2 * x) ** k for k, L in enumerate(coeff))


coeff = [L0, L1, L2]
assert isclose(ge, rk(1 - x, [L0, -L1, L2]), abs_tol=1e-10)
print(
    (
        f"Relabeling check: original="
        f"{ge:.9f}"
        f", wrong unchanged coefficients="
        f"{rk(1 - x, coeff):.9f}"
    )
)
# A symmetric even-order shape need not reduce to L0 alone.
print(
    (
        f"Even RK: gE(.25)="
        f"{rk(0.25, [-10000, 0, -10000]):.6f}"
        f", gE(.75)="
        f"{rk(0.75, [-10000, 0, -10000]):.6f}"
    )
)

xa, xb, xc = 0.2, 0.3, 0.5
binary = xa * xb * 1000 + xa * xc * 2000 + xb * xc * (-3000)
ternary = xa * xb * xc * (-100000)
assert isclose(binary, -190)
print(
    (
        f"Ternary (.2,.3,.5): binary sum="
        f"{binary:.6f}"
        f", correction="
        f"{ternary:.6f}"
        f", total="
        f"{binary + ternary:.6f}"
    )
)
for n in (4, 10):
    terms = [comb(n, k) / n**k for k in (2, 3, 4)]
    print(
        f"Equal coefficients, n={n}: total binary/ternary/quaternary weights="
        + ", ".join(f"{v:.9f}" for v in terms)
    )

# Pressure correction from a constant-parameter Murnaghan volume law.
V0 = 1e-5
B0 = 1e11
Bp = 4.0
dp = 1e9
v = V0 * (1 + Bp * dp / B0) ** (-1 / Bp)
gp = V0 * B0 / (Bp - 1) * ((1 + Bp * dp / B0) ** ((Bp - 1) / Bp) - 1)
print(
    (
        f"Pressure increase 1 GPa: gP="
        f"{gp:.9f}"
        f", V="
        f"{v:.12g}"
        f", incompressible V0 dp="
        f"{V0 * dp:.6f}"
    )
)

# Finite sampled phase candidates: exhaustive pairs plus a global LP dual
# certificate.
z = 0.35
A = 10000.0
points = [
    (phase, i / 20, A * (i / 20 - center) ** 2)
    for phase, center in [("alpha", 0.2), ("beta", 0.8)]
    for i in range(21)
]
candidates = []
for phase, xx, gg in points:
    if isclose(xx, z, abs_tol=1e-14):
        candidates.append((gg, (phase, xx, 1.0), (phase, xx, 0.0)))
for i, (p, u, gu) in enumerate(points):
    for q, vv, gv in points[i + 1:]:
        p, u, gu = points[i]
        if u > vv:
            p, u, gu, q, vv, gv = q, vv, gv, p, u, gu
        if u < vv and u <= z <= vv:
            f = (z - u) / (vv - u)
            candidates.append(
                ((1 - f) * gu + f * gv, (p, u, 1 - f), (q, vv, f))
            )
best = min(candidates, key=lambda t: t[0])
assert isclose(best[0], 0, abs_tol=1e-9)
assert all(gg >= -1e-10 for _, _, gg in points)
print(
    (
        f"Sampled equilibrium z=.35: g="
        f"{best[0]:.9f}"
        f", states="
        f"{best[1:]}"
        f"; dual line ell=0 certifies optimum"
    )
)


# Shared elemental reference shifts leave the optimizer unchanged.
def shift(xx):
    return 500 * (1 - xx) - 800 * xx


shifted_cost = best[1][2] * shift(best[1][1]) + best[2][2] * shift(best[2][1])
assert isclose(shifted_cost, shift(z), abs_tol=1e-10)
print(
    (
        f"Common elemental gauge: cost shift="
        f"{shifted_cost:.9f}"
        f", expected="
        f"{shift(z):.9f}"
        f"; phase amounts unchanged"
    )
)
