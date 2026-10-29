"""Original dimensionless teaching models; no assessed alloy database.

Run from this directory. Standard library only, no external effects.
Results test mathematical consistency, not real-material predictions.
"""

from math import isclose, log, sqrt
import json


def free_energy(x, interaction):
    """Molar mixing Gibbs energy divided by R*T for 0 <= x <= 1.

    x is component-B mole fraction; interaction is Omega/(R*T).
    Endmember reference energies are both zero in this toy model.
    """
    if not 0.0 <= x <= 1.0:
        raise ValueError("Mole fraction must lie in [0, 1].")
    if x in (0.0, 1.0):
        return 0.0
    return x * log(x) + (1.0 - x) * log(1.0 - x) + interaction * x * (1.0 - x)


def derivative(x, interaction):
    if not 0.0 < x < 1.0:
        raise ValueError("Derivative requires an interior composition.")
    return log(x / (1.0 - x)) + interaction * (1.0 - 2.0 * x)


def coexistence(interaction):
    """Symmetric regular-solution binodal endpoints for interaction > 2.

    Equal endpoint energy makes the common-tangent slope zero; solve
    g'(x)=0 in (0, lower spinodal), not the unstable central root 1/2.
    """
    if interaction <= 2.0:
        raise ValueError("Distinct coexistence requires Omega/(R*T) > 2.")
    spinodal = (1.0 - sqrt(1.0 - 2.0 / interaction)) / 2.0
    lo, hi = 1.0e-15, spinodal
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if derivative(mid, interaction) > 0.0:
            hi = mid
        else:
            lo = mid
    x_alpha = (lo + hi) / 2.0
    return x_alpha, 1.0 - x_alpha


def phase_fractions(overall, alpha, beta):
    """Mole fractions of phases on a two-phase tie line, from B balance."""
    if not 0.0 <= alpha < beta <= 1.0 or not alpha <= overall <= beta:
        raise ValueError("Ordered endpoints must enclose overall composition.")
    f_beta = (overall - alpha) / (beta - alpha)
    return 1.0 - f_beta, f_beta


def check_model():
    interaction = 3.0
    alpha, beta = coexistence(interaction)
    spinodal = (1.0 - sqrt(1.0 - 2.0 / interaction)) / 2.0
    assert 0.0 < alpha < spinodal < 0.5 < 1.0 - spinodal < beta < 1.0
    assert abs(derivative(alpha, interaction)) < 1.0e-12
    assert abs(derivative(beta, interaction)) < 1.0e-12
    energy = free_energy(alpha, interaction)
    assert isclose(energy, free_energy(beta, interaction), abs_tol=1.0e-12)
    residual = min(
        free_energy(i / 5000.0, interaction) - energy for i in range(5001)
    )
    assert residual >= -1.0e-12
    checks = 0
    for i in range(101):
        overall = alpha + (beta - alpha) * i / 100.0
        fa, fb = phase_fractions(overall, alpha, beta)
        assert isclose(fa + fb, 1.0, abs_tol=1.0e-12)
        assert isclose(fa * alpha + fb * beta, overall, abs_tol=1.0e-12)
        checks += 1
    return {
        "model": "symmetric regular solution; synthetic teaching model",
        "interaction_Omega_over_RT": interaction,
        "binodal": [alpha, beta],
        "spinodal": [spinodal, 1.0 - spinodal],
        "common_tangent_slope": 0.0,
        "coexistence_g_over_RT": energy,
        "mass_balance_checks": checks,
        "sampled_supporting_line_min_residual": residual,
        "phase_mole_fractions_at_x_0_4": phase_fractions(0.4, alpha, beta),
        "limitations": (
            "Finite grid check supplements model derivation; "
            "not a measured alloy."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(check_model(), indent=2))
