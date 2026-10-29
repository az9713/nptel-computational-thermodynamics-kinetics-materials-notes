"""Declared ideal models; standard library only; no experimental data."""

from math import exp, log, factorial, isclose, sqrt

R = 8.31446261815324
C = 1.5 * R
TREF, VREF = 300.0, 0.025


def entropy(t, v):
    return C * log(t / TREF) + R * log(v / VREF)


def energy(s, v):
    return C * TREF * exp(s / C) * (VREF / v) ** (R / C)


def availability(t, v, t0, p0):
    return C * t + p0 * v - t0 * entropy(t, v)


def phase_difference(t):
    dh = 10000 + 20 * (t - 1000)
    ds = 10 + 20 * log(t / 1000)
    return dh, ds, dh - t * ds


def main():
    wi, wf = 56 * 4, 20 * 20
    hot, cold = log(20 / 56), log(20 / 4)
    assert isclose(hot + cold, log(wf / wi))
    print(
        (
            f"Four-mode transfer: dS_hot/k="
            f"{hot:.9f}"
            f"; cold="
            f"{cold:.9f}"
            f"; total="
            f"{hot + cold:.9f}"
        )
    )
    # Correlated equiprobable pairs00 and11: two fair marginals but only two
    # joint states.
    joint = log(2)
    marginal_sum = 2 * log(2)
    mutual = marginal_sum - joint
    print(
        (
            f"Correlated binary pair: S_joint/k="
            f"{joint:.9f}"
            f"; sum_marginals="
            f"{marginal_sum:.9f}"
            f"; mutual="
            f"{mutual:.9f}"
        )
    )
    # The exact one-mole gas fundamental relation has positive d2U/dS2=T/C.
    t, v = 450.0, 0.04
    s = entropy(t, v)
    h = 1e-4
    t_num = (energy(s + h, v) - energy(s - h, v)) / (2 * h)
    p_num = -(energy(s, v + 1e-7) - energy(s, v - 1e-7)) / (2e-7)
    assert isclose(t_num, t, rel_tol=1e-8)
    assert isclose(p_num, R * t / v, rel_tol=1e-8)
    a = energy(s, v) - t * s
    assert isclose(a, C * t - t * s, abs_tol=1e-10)
    print(
        (
            f"Ideal fundamental relation: S="
            f"{s:.9f}"
            f" J/K; T="
            f"{t_num:.6f}"
            f" K; p="
            f"{p_num:.6f}"
            f" Pa; A="
            f"{a:.6f}"
            f" J"
        )
    )
    # Availability equals G at endpoints whose T,p match the environment.
    t0, p0 = 300.0, 100000.0
    vi = R * 600 / p0
    vf = R * t0 / p0
    delta_phi = availability(t0, vf, t0, p0) - availability(600, vi, t0, p0)
    # Constant external pressure/no useful work: Q=dH=Cp*(Tf-Ti).
    dh = (C + R) * (300 - 600)
    ds = (C + R) * log(300 / 600)
    gen = ds - dh / t0
    assert isclose(delta_phi, -t0 * gen, abs_tol=1e-9)
    print(
        (
            f"600-to300 K cooling at100 kPa: Q="
            f"{dh:.6f}"
            f" J; dS_gas="
            f"{ds:.9f}"
            f" J/K; generation="
            f"{gen:.9f}"
            f" J/K; dPhi="
            f"{delta_phi:.6f}"
            f" J"
        )
    )
    print(
        "Temperature-dependent phase differences: "
        + "; ".join(
            (
                f"T="
                f"{t}"
                f": dH="
                f"{phase_difference(t)[0]:.6f}"
                f", dS="
                f"{phase_difference(t)[1]:.9f}"
                f", dG="
                f"{phase_difference(t)[2]:.6f}"
            )
            for t in (500, 1000, 1500)
        )
    )
    # One quantum oscillator, excitation energies j*k*theta, excluding
    # zero-point shift.
    x = 2.0
    z = 1 / (1 - exp(-x))
    nbar = 1 / (exp(x) - 1)
    s_mode = x * nbar + log(z)
    assert isclose(s_mode, x / (exp(x) - 1) - log(1 - exp(-x)))
    print(
        (
            f"Oscillator theta/T=2: Z="
            f"{z:.9f}"
            f"; mean_quanta="
            f"{nbar:.9f}"
            f"; S/k="
            f"{s_mode:.9f}"
        )
    )
    # Analytic entropy extraction: small noise in G can have amplified slope.
    # G(T)=a+b*T-c*T*log(T/Tref), so S=-b+c*(log(T/Tref)+1).
    b, c, t = -10.0, 20.0, 600.0

    def g(temp):
        return 1000 + b * temp - c * temp * log(temp / 300)

    ds_exact = -b + c * (log(t / 300) + 1)
    eps = 0.1
    ds_num = -(g(t + eps) - g(t - eps)) / (2 * eps)
    assert isclose(ds_num, ds_exact, rel_tol=1e-8)
    slope_noise = sqrt(2) * 0.5 / (2 * eps)
    print(
        (
            f"Analytic entropy slope: exact="
            f"{ds_exact:.9f}"
            f", finite_difference="
            f"{ds_num:.9f}"
            f" J/(mol K); sigma_G=0.5 J/mol gives sigma_slope="
            f"{slope_noise:.9f}"
        )
    )
    assert phase_difference(500)[2] > 0 and phase_difference(1500)[2] < 0
    print(
        (
            "All fundamental-relation, availability, counting and finite-diff"
            "erence checks passed."
        )
    )


if __name__ == "__main__":
    main()
