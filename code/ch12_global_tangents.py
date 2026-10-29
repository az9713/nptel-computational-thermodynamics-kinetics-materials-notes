"""Exact declared parabolic certificates and bounded regular-solution c\
hecks."""

from math import log, sqrt, isclose

R = 8.31446261815324
A, LEFT, RIGHT = 10000.0, 0.2, 0.8


def ga(x):
    return A * (x - LEFT) ** 2


def gb(x, delta=0):
    return A * (x - RIGHT) ** 2 + delta


def pair(delta=0):
    slope = delta / (RIGHT - LEFT)
    xa = LEFT + slope / (2 * A)
    xb = RIGHT + slope / (2 * A)
    intercept = ga(xa) - slope * xa
    return xa, xb, slope, intercept


def g_regular(x, t, omega=12000):
    return R * t * (x * log(x) + (1 - x) * log(1 - x)) + omega * x * (1 - x)


def bisect(fn, a, b):
    fa = fn(a)
    fb = fn(b)
    assert fa * fb < 0
    for _ in range(100):
        c = (a + b) / 2
        fc = fn(c)
        if fa * fc <= 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return (a + b) / 2


def main():
    for delta in (0, 500):
        xa, xb, m, c = pair(delta)
        assert 0 < xa < xb < 1
        for j in range(1001):
            x = j / 1000
            # Both exact residual identities are global certificates, not mere
            # grid claims.
            assert isclose(
                ga(x) - (c + m * x), A * (x - xa) ** 2, abs_tol=1e-9
            )
            assert isclose(
                gb(x, delta) - (c + m * x), A * (x - xb) ** 2, abs_tol=1e-9
            )
        xt0 = 0.5 + delta / (2 * A * (RIGHT - LEFT))
        assert isclose(ga(xt0), gb(xt0, delta), abs_tol=1e-9)
        print(
            (
                f"Parabolic delta="
                f"{delta}"
                f": coexistence="
                f"{xa:.9f}"
                f","
                f"{xb:.9f}"
                f"; slope="
                f"{m:.9f}"
                f"; T0 composition="
                f"{xt0:.9f}"
            )
        )
    z = 0.3
    fraction = (z - LEFT) / (RIGHT - LEFT)
    assert isclose((1 - fraction) * LEFT + fraction * RIGHT, z)
    print(
        (
            f"At overallxB0.3: beta mole fraction="
            f"{fraction:.9f}"
            f"; finite release="
            f"{ga(z):.9f}"
            f" J/mol_initial"
        )
    )

    # Old alpha/beta equality alone misses a new competing phase.
    def gamma(x):
        return A * (x - 0.5) ** 2 - 200

    assert gamma(0.5) < 0
    print(
        (
            f"Third-phase counterexample: old common line=0; competing gam"
            f"ma at0.5="
            f"{gamma(0.5):.6f}"
            f" J/mol"
        )
    )
    parent_slope = 2 * A * (z - LEFT)

    def tangent(y):
        return ga(z) + parent_slope * (y - z)

    ystar = RIGHT + parent_slope / (2 * A)
    drive = tangent(ystar) - gb(ystar)

    def infant(eta):
        xp = (z - eta * ystar) / (1 - eta)
        return (1 - eta) * ga(xp) + eta * gb(ystar)

    eta = 1e-7
    derivative = (infant(eta) - ga(z)) / eta
    assert isclose(derivative, -drive, abs_tol=0.001)
    print(
        (
            f"Infinitesimal product: xB="
            f"{ystar:.9f}"
            f"; parent-tangent release="
            f"{drive:.9f}"
            f" J/mol_product; finite_difference_derivative="
            f"{derivative:.6f}"
        )
    )
    # A separate full-progress path keeps productcomposition at its
    # final0.8value.
    d = z - LEFT
    h = RIGHT - LEFT

    def progress_g(eta):
        return A * (d - h * eta) ** 2 / (1 - eta)

    def affinity(eta):
        return A * (
            2 * h * (d - h * eta) / (1 - eta)
            - (d - h * eta) ** 2 / (1 - eta) ** 2
        )

    n = 1000
    step = fraction / n
    integral = (
        step
        * (
            affinity(0)
            + affinity(fraction)
            + sum(
                (4 if j % 2 else 2) * affinity(j * step) for j in range(1, n)
            )
        )
        / 3
    )
    assert isclose(
        integral, progress_g(0) - progress_g(fraction), abs_tol=1e-9
    )
    print(
        (
            f"Full constrained path: initial affinity="
            f"{affinity(0):.9f}"
            f" J/mol_product; integral="
            f"{integral:.9f}"
            f" J/mol_initial; final affinity="
            f"{affinity(fraction):.9f}"
        )
    )
    gc, xc = -100.0, 0.5
    # Symmetric tangent from alpha to the line compound.
    xa = xc - sqrt((xc - LEFT) ** 2 - gc / A)
    ml = 2 * A * (xa - LEFT)
    assert isclose(ga(xa) + ml * (xc - xa), gc, abs_tol=1e-9)
    print(
        (
            f"Line compound: x_alpha="
            f"{xa:.9f}"
            f"; x_beta="
            f"{1 - xa:.9f}"
            f"; admissible purecompound support slopes=["
            f"{ml:.9f}"
            f","
            f"{-ml:.9f}"
            f"]"
        )
    )
    t, omega = 600.0, 12000.0
    tc = omega / (2 * R)
    spin = (1 - sqrt(1 - 2 * R * t / omega)) / 2
    binodal = bisect(
        lambda x: R * t * log(x / (1 - x)) + omega * (1 - 2 * x),
        1e-12,
        0.499999,
    )
    assert binodal < spin < 0.5
    print(
        (
            f"Regular model: criticalT="
            f"{tc:.9f}"
            f" K; binodal="
            f"{binodal:.9f}"
            f","
            f"{1 - binodal:.9f}"
            f"; spinodal="
            f"{spin:.9f}"
            f","
            f"{1 - spin:.9f}"
        )
    )
    # Check a convex local branch that is still above the global two-phase
    # line.
    test = 0.5 * (binodal + spin)
    curvature = R * t / (test * (1 - test)) - 2 * omega
    gap = g_regular(test, t) - g_regular(binodal, t)
    assert curvature > 0 and gap > 0
    print(
        (
            f"Metastable local test: xB="
            f"{test:.9f}"
            f"; curvature="
            f"{curvature:.9f}"
            f"; above equilibrium line="
            f"{gap:.9f}"
            f" J/mol"
        )
    )
    print(
        (
            "Generic intensive counts: binary2phase F=2, fixedp F=1; binary3p"
            "hase fixedp F=0; ternary2phase fixedTp F=1"
        )
    )
    print(
        (
            "All exact supporting-line, conservation, affinity and bounded ro"
            "ot checks passed."
        )
    )


if __name__ == "__main__":
    main()
