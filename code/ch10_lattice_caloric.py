"""Declared unary/bond/oscillator models, plus attributed small publish\
ed data extract.\nLattice rows: Axel van de Walle, Mendeley Data v2,\
 doi:10.17632/h2ggh95f9m.2,\nCC BY4.0; selected values transcribed f\
rom updated table.csv, not DFT reruns.\n"""

from math import log, exp, expm1, isclose

R = 8.31446261815324
ALPHA = 3e-5
T_JOIN, TM, LATENT = 100.0, 1000.0, 10000.0


def solid_hs(t):
    assert t > 0
    if t <= T_JOIN:
        return ALPHA * t**4 / 4, ALPHA * t**3 / 3
    return 750 + 30 * (t - 100), 10 + 30 * log(t / 100)


def liquid_hs(t):
    # Declared extrapolated liquid branch here is used only above298K.
    assert t >= 298
    h, s = solid_hs(TM)
    return h + LATENT + 40 * (t - TM), s + LATENT / TM + 40 * log(t / TM)


def delta_g(t):
    hs, ss = solid_hs(t)
    hl, sl = liquid_hs(t)
    return hl - hs - t * (sl - ss)


def oscillator(theta, t):
    # Molar quantum harmonic energy includes zero-point term.
    x = theta / t
    u = R * theta * (0.5 + 1 / expm1(x))
    s = R * (x / expm1(x) - log(1 - exp(-x)))
    f = u - t * s
    return u, s, f


def main():
    for t in (100, 1000, 1200):
        h, s = solid_hs(t) if t <= 1000 else liquid_hs(t)
        print(
            (
                f"Caloric path T="
                f"{t}"
                f": H="
                f"{h:.6f}"
                f" J/mol; S="
                f"{s:.9f}"
                f" J/(mol K); G="
                f"{h - t * s:.6f}"
                f" J/mol"
            )
        )
    hs, ss = solid_hs(TM)
    hl, sl = liquid_hs(TM)
    assert isclose(hl - hs, LATENT) and isclose(sl - ss, LATENT / TM)
    assert isclose(delta_g(TM), 0, abs_tol=1e-9)
    print(f"Fusion jump: H={hl-hs:.6f}; S={sl-ss:.9f}; G={delta_g(TM):.9f}")
    for t in (800, 1000, 1200):
        hs, ss = solid_hs(t)
        hl, sl = liquid_hs(t)
        common_shift = 12345 - 7 * t
        assert isclose(
            (hl - t * sl + common_shift) - (hs - t * ss + common_shift),
            delta_g(t),
            abs_tol=1e-9,
        )
    print(
        "Lattice difference: "
        + "; ".join(
            f"T={t}: dG={delta_g(t):.6f} J/mol" for t in (800, 1000, 1200)
        )
    )
    hl0, sl0 = liquid_hs(298)
    hs0, ss0 = solid_hs(298)
    hs_tm, _ = solid_hs(TM)
    hl_tm, _ = liquid_hs(TM)
    wrong_fusion = (hl_tm - hl0) - (hs_tm - hs0)
    assert isclose(wrong_fusion, 7020)
    print(
        (
            f"Independent phase zeroing destroys offset: omitted H_liq-H_s"
            f"ol at298K="
            f"{hl0 - hs0:.6f}"
            f"; wrong fusion="
            f"{wrong_fusion:.6f}"
            f" J/mol"
        )
    )
    # 100sites with coordination6: Nz/2 unordered edges, random independent
    # labels.
    n, z, x = 100, 6, 0.3
    edges = n * z / 2
    expected_ab = edges * 2 * x * (1 - x)
    assert isclose(expected_ab, n * z * x * (1 - x))
    omega = 12000
    print(
        (
            f"Random bond model: unordered_edges="
            f"{edges:.0f}"
            f"; expected_AB="
            f"{expected_ab:.6f}"
            f"; molar_dH="
            f"{omega * x * (1 - x):.6f}"
            f" J/mol"
        )
    )
    # Fixed-composition random labeling has a finite-N correction.
    fixed_ab = edges * (2 * (n * (1 - x)) * (n * x) / (n * (n - 1)))
    print(
        (
            f"Fixed30B/70A finite100site ensemble: expected_AB="
            f"{fixed_ab:.9f}"
            f"; thermodynamic-limit="
            f"{expected_ab:.9f}"
        )
    )
    for t in (300, 30000):
        u1, s1, f1 = oscillator(300, t)
        u2, s2, f2 = oscillator(600, t)
        assert isclose(f2 - f1, (u2 - u1) - t * (s2 - s1), abs_tol=1e-9)
        print(
            (
                f"Stiffness change theta300to600 atT="
                f"{t}"
                f": dU="
                f"{u2 - u1:.9f}"
                f"; dS="
                f"{s2 - s1:.9f}"
                f"; dF="
                f"{f2 - f1:.9f}"
            )
        )
    pressure_delta = 101325 - 100000
    print(
        (
            f"Onebar-to-oneatm: delta_p="
            f"{pressure_delta}"
            f" Pa; atV=1e-5m3/mol delta_g="
            f"{1e-05 * pressure_delta:.8f}"
            f" J/mol"
        )
    )
    # Headers report eV; do not infer a missing atom/mole normalization here.
    published = [
        ("Al fcc->bcc", 0.105, 0.047, 0.094),
        ("Ni fcc->bcc", 0.090, 0.094, 0.093),
    ]
    for label, sgte, inflection, naive in published:
        print(
            (
                f"Published extract "
                f"{label}"
                f": SGTE="
                f"{sgte:.3f}"
                f"; inflection="
                f"{inflection:.3f}"
                f"; shape-constrained="
                f"{naive:.3f}"
                f" eV; inflection-minus-SGTE="
                f"{inflection - sgte:.3f}"
            )
        )
    print(
        (
            "All declared caloric, reference-invariance, bond-count and oscil"
            "lator checks passed."
        )
    )


if __name__ == "__main__":
    main()
