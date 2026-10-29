"""Synthetic property assessments and strength checks; no real alloy fit.
Run with Python 3.13. Inputs are stated teaching parameters in SI units,
except the explicitly labelled MPa/nm strength model. No files are written.
"""

from math import exp, isclose, log, sqrt

GAS_R = 8.31446261815324
T0, P0 = 1000.0, 1e5
V0, ALPHA0, KAPPA0, CP = 1e-5, 6e-5, 1e-11, 25.0


def gibbs(temperature, pressure):
    if temperature <= 0:
        raise ValueError("Positive temperature required")
    q = pressure - P0
    thermal = -CP * (temperature * log(temperature / T0) - temperature + T0)
    return thermal + V0 * (
        (1 + ALPHA0 * (temperature - T0)) * q - KAPPA0 * q * q / 2
    )


def response(temperature, pressure):
    if temperature <= 0:
        raise ValueError("Positive temperature required")
    factor = 1 + ALPHA0 * (temperature - T0) - KAPPA0 * (pressure - P0)
    if factor <= 0:
        raise ValueError("Outside positive-volume domain of the local model")
    volume = V0 * factor
    alpha, kappa = ALPHA0 / factor, KAPPA0 / factor
    cv = CP - temperature * volume * alpha * alpha / kappa
    if cv <= 0:
        raise ValueError("Outside thermal-stability domain of the local model")
    return volume, alpha, kappa, 1 / kappa, cv


volume, alpha, kappa, modulus, cv = response(1100.0, P0 + 1e9)
hp, ht = 1e6, 0.1
vp_fd = (gibbs(1100.0, P0 + 1e9 + hp) - gibbs(1100.0, P0 + 1e9 - hp)) / (
    2 * hp
)
gpp_fd = (
    gibbs(1100.0, P0 + 1e9 + hp)
    - 2 * gibbs(1100.0, P0 + 1e9)
    + gibbs(1100.0, P0 + 1e9 - hp)
) / (hp * hp)
cp_fd = (
    -1100
    * (
        gibbs(1100 + ht, P0 + 1e9)
        - 2 * gibbs(1100.0, P0 + 1e9)
        + gibbs(1100 - ht, P0 + 1e9)
    )
    / (ht * ht)
)
assert isclose(vp_fd, volume, rel_tol=1e-9)
assert isclose(-gpp_fd / volume, kappa, rel_tol=1e-6)
assert isclose(cp_fd, CP, rel_tol=1e-5)
print(f"Local molar volume at 1100 K, pressure increment 1 GPa: {volume:.9e}")
print(f"Volumetric expansivity and compressibility: {alpha:.9e}, {kappa:.9e}")
print(
    f"Isothermal bulk modulus and Cv: {modulus/1e9:.6f} GPa, "
    f"{cv:.6f} J/(mol K)"
)
print(f"Finite-difference volume and Cp: {vp_fd:.9e}, {cp_fd:.6f}")
print(f"Maxwell entropy-pressure derivative: {-V0*ALPHA0:.9e} m3/(mol K)")

# Independent log-viscosity model: x is B atomic fraction in the liquid.
ETA_A, ETA_B = 0.004, 0.002  # Pa s at T0
E_A, E_B = 20000.0, 15000.0  # J/mol
L0, L1 = 6000.0, -2000.0  # J/mol; not Gibbs interaction parameters


def viscosity(x, temperature=T0, l0=L0, l1=L1):
    if not 0 <= x <= 1 or temperature <= 0:
        raise ValueError("Require composition in [0,1] and temperature > 0")
    ya = log(ETA_A) + E_A / GAS_R * (1 / temperature - 1 / T0)
    yb = log(ETA_B) + E_B / GAS_R * (1 / temperature - 1 / T0)
    interaction = x * (1 - x) * (l0 + l1 * (1 - 2 * x)) / (GAS_R * temperature)
    return exp((1 - x) * ya + x * yb + interaction)


def reduced_interaction(x):
    pure = (1 - x) * log(ETA_A) + x * log(ETA_B)
    return GAS_R * T0 * (log(viscosity(x)) - pure) / (x * (1 - x))


eta_mid = viscosity(0.5)
assert isclose(viscosity(0), ETA_A, rel_tol=1e-14)
assert isclose(viscosity(1), ETA_B, rel_tol=1e-14)
assert isclose(viscosity(0.5, l1=1e6), eta_mid, rel_tol=1e-14)
d25, d75 = reduced_interaction(0.25), reduced_interaction(0.75)
fit0, fit1 = (d25 + d75) / 2, d25 - d75
assert isclose(fit0, L0, abs_tol=1e-8)
assert isclose(fit1, L1, abs_tol=1e-8)
rho_ref = (0.027 + 0.045) / 2 / V0  # Reference equimolar mass/volume: kg/m3
print(
    f"Viscosity at x=.25, .5, .75: {viscosity(.25):.9f}, "
    f"{eta_mid:.9f}, {viscosity(.75):.9f} Pa s"
)
print(
    f"Equimolar no-excess viscosity and multiplier: "
    f"{sqrt(ETA_A*ETA_B):.9f}, {eta_mid/sqrt(ETA_A*ETA_B):.9f}"
)
print(f"Two-composition reconstructed L0, L1: {fit0:.6f}, {fit1:.6f} J/mol")
print(
    f"Reference density and kinematic viscosity: "
    f"{rho_ref:.6f} kg/m3, {eta_mid/rho_ref:.9e} m2/s"
)
# Same local conductivities, different imposed one-dimensional geometry.
k_parallel, k_series = 0.5 * 100 + 0.5 * 10, 1 / (0.5 / 100 + 0.5 / 10)
print(
    f"Parallel/series thermal conductivity: "
    f"{k_parallel:.6f}, {k_series:.6f} W/(m K)"
)

# (Fe,Ti)1(C,Va)1: hypothetical endmember strengths, not measured values.
y_ti, y_c = 0.2, 0.1
weights = (
    (1 - y_ti) * (1 - y_c),
    y_ti * (1 - y_c),
    (1 - y_ti) * y_c,
    y_ti * y_c,
)
endmember = (50.0, 100.0, 200.0, 500.0)  # Fe:Va, Ti:Va, Fe:C, Ti:C in MPa
intrinsic = sum(w * s for w, s in zip(weights, endmember))
assert isclose(sum(weights), 1.0, abs_tol=1e-15)
x_c, x_ti = y_c / (1 + y_c), y_ti / (1 + y_c)
print("Endmember weights: " + ", ".join(f"{w:.6f}" for w in weights))
print(
    f"Intrinsic strength and atomic C/Ti fractions: "
    f"{intrinsic:.6f} MPa, {x_c:.9f}, {x_ti:.9f}"
)
# These additional constitutive choices are independent teaching parameters.
grain_size, hp_coefficient = 25e-6, 0.08  # m and MPa sqrt(m)
matrix_solute, ss_coefficient = 0.02, 400.0  # atomic fraction and MPa
gb_strength = hp_coefficient / sqrt(grain_size)
ss_strength = ss_coefficient * matrix_solute ** (2 / 3)
base_strength = intrinsic + gb_strength + ss_strength
print(
    f"Grain and matrix-solute increments: "
    f"{gb_strength:.6f}, {ss_strength:.6f} MPa"
)
print(f"Linear baseline without precipitation: {base_strength:.6f} MPa")
print(
    f"Using bulk solute .04 instead of matrix .02 gives: "
    f"{ss_coefficient*.04**(2/3):.6f} MPa"
)

# Stipulated constant-volume-fraction cutting/bypass envelope, radius in nm.
CUT, BYPASS = 40.0, 800.0  # Both coefficients in MPa; u = r/(1 nm)


def increments(radius_nm):
    if radius_nm <= 0:
        raise ValueError("Positive radius required")
    return CUT * sqrt(radius_nm), BYPASS / radius_nm


critical = (BYPASS / CUT) ** (2 / 3)
peak = CUT * sqrt(critical)
assert isclose(*increments(critical), rel_tol=1e-14)
assert min(increments(critical / 2)) < peak
assert min(increments(2 * critical)) < peak
print(f"Toy crossover radius and increment: {critical:.9f} nm, {peak:.6f} MPa")
print(f"Toy total strength at the crossover: {base_strength+peak:.6f} MPa")
local = [increments(r) for r in (2.0, 20.0)]
average_minimum = sum(min(pair) for pair in local) / 2
minimum_average = min(sum(pair[i] for pair in local) / 2 for i in (0, 1))
mean_radius_value = min(increments(11.0))
assert average_minimum <= minimum_average
print(
    f"Distribution operation checks, average(min), min(average), "
    f"min(at mean r): {average_minimum:.6f}, {minimum_average:.6f}, "
    f"{mean_radius_value:.6f} MPa"
)
print(
    f"Two stipulated 100 MPa increments: linear={200.:.6f}, "
    f"quadrature={sqrt(100**2+100**2):.6f} MPa"
)
