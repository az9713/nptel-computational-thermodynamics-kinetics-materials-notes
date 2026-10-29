"""Original size, transport, electrochemical and database-scope audits.
Synthetic parameters only: not an assessed alloy/corrosion/process simulation.
"""

from math import pi, exp, log, log10, comb, isclose

R, F = 8.31446261815324, 96485.33212331002
vm, gamma, radius = 1e-5, 1.0, 10e-9
n = 4 * pi * radius**3 / (3 * vm)


def area(moles):
    return 4 * pi * (3 * vm * moles / (4 * pi)) ** (2 / 3)


average = gamma * area(n) / n
partial = 2 * gamma * vm / radius
step = n * 1e-5
finite = gamma * (area(n + step) - area(n - step)) / (2 * step)
assert isclose(average, 3000.0, rel_tol=1e-12)
assert isclose(partial, 2000.0, rel_tol=1e-12)
assert isclose(finite, partial, rel_tol=1e-9)
print(
    (
        f"Spherical surface contribution, total/n and derivative: "
        f"{average:.6f}"
        f", "
        f"{partial:.6f}"
        f" J/mol"
    )
)
print(f"Centered finite-difference derivative: {finite:.6f} J/mol")
# Equal volumes and equal constrained radius, constant surface energies.
for rnm in (5.0, 10.0, 50.0):
    shift = 3 * vm * (0.6 - 1.0) / (rnm * 1e-9)
    tm = 1000.0 + shift / 10.0
    assert abs(10 * (1000 - tm) + shift) < 1e-9
    print(f"Fixed-radius two-state crossing, radius {rnm:g} nm: {tm:.6f} K")

# Toy parallel/series scalar paths; geometric weights are stipulated.
T, Q = 1000.0, 200000.0
fgb, fp = 0.5e-9 / (50e-9), 1e15 * (0.25e-9) ** 2
fb = 1 - fgb - fp
assert min(fgb, fp, fb) >= 0
mgb, mp = exp((1 - 0.5) * Q / (R * T)), exp((1 - 0.8) * Q / (R * T))
parallel = fb + fgb * mgb + fp * mp
series = 1 / (fb + fgb / mgb + fp / mp)
assert series <= parallel and fb + fgb + fp == 1.0
print(f"Toy path fractions bulk/GB/pipe: {fb:.8f}, {fgb:.8f}, {fp:.8f}")
print(
    (
        f"Toy effective/bulk mobility, parallel and series: "
        f"{parallel:.6f}"
        f", "
        f"{series:.6f}"
    )
)

# Reduction Ox + mH+ + ne- = Red; all other activities unity.
T = 298.15
slope = R * T * log(10) / F
print(f"Nernst decade factor at 298.15 K: {slope:.9f} V")
for ph in (0.0, 7.0, 14.0):
    potential = 0.8 - slope * ph  # m=n=2
    print(
        f"Illustrative m/n=1 reduction boundary, pH {ph:g}: {potential:.9f} V"
    )
metal = -0.2 + slope / 2 * log10(1e-6)
print(f"Illustrative M2+/M boundary at ion activity 1e-6: {metal:.9f} V")

# M + .5 O2 -> MO, unit condensed activities and ideal oxygen gas.
T, dg0 = 1500.0, -200000.0
logp = 2 * dg0 / (R * T * log(10))
assert abs(dg0 - 0.5 * R * T * log(10**logp)) < 1e-8
print(f"Toy oxide equilibrium log10(pO2/pstandard): {logp:.9f}")
print(
    "Five-component binary/ternary subsystem counts: "
    f"{comb(5, 2)}, {comb(5, 3)}"
)

# Two ternary models agree on every binary edge but differ in the interior.
for sign in (-1.0, 1.0):
    interaction = sign * 27000.0 * (1 / 3) ** 3
    assert sign * 27000.0 * 0.0 * 0.4 * 0.6 == 0.0
    print(
        (
            f"Edge-invisible ternary interaction at center: "
            f"{interaction:.6f}"
            f" J/mol"
        )
    )
# Closed binary: pure A/B at zero, an optional AB line phase at x=.5.
z = 0.4
old = 0.0
sigma_fraction = 2 * min(z, 1 - z)
new = -500 * sigma_fraction
assert new <= old
print(
    (
        f"Missing-phase toy at x=.4, restricted/enlarged envelope: "
        f"{old:.6f}"
        f", "
        f"{new:.6f}"
        f" J/mol"
    )
)
# Stipulated bounded model discrepancy; not inferred confidence levels.
margin, error = -300.0, 500.0
print(
    (
        f"Energy-difference interval under stipulated bound: ["
        f"{margin - error:.1f}"
        f", "
        f"{margin + error:.1f}"
        f"] J/mol"
    )
)
