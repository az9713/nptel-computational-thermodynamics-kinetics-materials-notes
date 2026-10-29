"""Declared monatomic gas paths and an independent harmonic-quench integral."""

from math import log, exp, sqrt, pi, isclose

R = 8.31446261815324
N, T1, P1 = 1.0, 300.0, 100000.0
CV, CP = 1.5 * R, 2.5 * R
GAMMA = CP / CV
V1 = N * R * T1 / P1


def stepped_isothermal(steps, ratio=2.0):
    # Increasing equilibrium volumes; each expansion acts against its final
    # load.
    volumes = [V1 * (1 + (ratio - 1) * i / steps) for i in range(steps + 1)]
    return sum(
        N * R * T1 / volumes[i] * (volumes[i] - volumes[i - 1])
        for i in range(1, steps + 1)
    )


def simpson(f, a, b, n=6000):
    step = (b - a) / n
    return (
        step
        / 3
        * (
            f(a)
            + f(b)
            + 4 * sum(f(a + i * step) for i in range(1, n, 2))
            + 2 * sum(f(a + i * step) for i in range(2, n, 2))
        )
    )


def main():
    rev = N * R * T1 * log(2)
    one = N * R * T1 / 2
    compression = -N * R * T1
    assert (
        one
        < stepped_isothermal(2)
        < stepped_isothermal(10)
        < stepped_isothermal(100)
        < rev
    )
    beta = 0.5
    trev = T1 * beta ** ((GAMMA - 1) / GAMMA)
    vrev = N * R * trev / (beta * P1)
    wrev = N * CV * (T1 - trev)
    tirr = T1 * (CV + beta * R) / CP
    virr = N * R * tirr / (beta * P1)
    wirr = beta * P1 * (virr - V1)
    assert isclose(N * CV * (tirr - T1), -wirr, abs_tol=1e-10)
    assert isclose(P1 * V1**GAMMA, beta * P1 * vrev**GAMMA, rel_tol=1e-12)
    assert tirr > trev and virr > vrev and wirr < wrev
    # A subsequent one-step adiabatic compression against the original load.
    treturn = (CV * tirr + P1 * virr / N) / CP
    vreturn = N * R * treturn / P1
    wreturn = P1 * (vreturn - virr)
    assert isclose(N * CV * (treturn - tirr), -wreturn, abs_tol=1e-10)
    assert treturn > T1 and vreturn > V1
    entropy_irr = N * CV * log(tirr / T1) + N * R * log(virr / V1)
    assert entropy_irr > 0

    # x scaled by initial Gaussian width; trap stiffness suddenly multiplies by
    # four.
    def gaussian(q):
        return exp(-q * q / 2) / sqrt(2 * pi)

    mean_work = simpson(lambda q: 1.5 * q * q * gaussian(q), -10, 10)
    exponential = simpson(lambda q: exp(-1.5 * q * q) * gaussian(q), -10, 10)
    assert isclose(mean_work, 1.5, abs_tol=1e-11)
    assert isclose(exponential, 0.5, abs_tol=1e-11)
    assert mean_work > -log(exponential)
    print(
        "Isothermal doubling: reversible %.6f J; one-step %.6f J." % (rev, one)
    )
    print(
        "Two / ten / hundred load steps: "
        + ", ".join("%.6f" % stepped_isothermal(k) for k in [2, 10, 100])
        + " J."
    )
    print(
        (
            "One-step isothermal compression W_by = %.6f J; full return cycle"
            " W_by = %.6f J."
        )
        % (compression, one + compression)
    )
    print(
        "Adiabatic to 50 kPa: reversible T %.6f K, V/V1 %.6f, W_by %.6f J."
        % (trev, vrev / V1, wrev)
    )
    print(
        "Adiabatic to 50 kPa: one-step T %.6f K, V/V1 %.6f, W_by %.6f J."
        % (tirr, virr / V1, wirr)
    )
    print(
        "Adiabatic one-step return to 100 kPa: T %.6f K; V/V1 %.6f."
        % (treturn, vreturn / V1)
    )
    print("Adiabatic one-step gas entropy change = %.9f J/K." % entropy_irr)
    print(
        (
            "Harmonic stiffness x4 quench: mean W_on/kBT %.9f; exp average %."
            "9f; Delta F/kBT %.9f."
        )
        % (mean_work, exponential, -log(exponential))
    )
    print(
        (
            "Checks: common isothermal endpoints, monotone load sums, both ad"
            "iabatic balances and quench integrals PASS."
        )
    )


if __name__ == "__main__":
    main()
