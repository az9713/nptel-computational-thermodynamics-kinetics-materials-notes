"""Original CEF, association and weighted-assessment examples; no asses\
sed alloy."""

from math import log, sqrt, isclose, comb

R = 8.31446261815324


def q(y):
    return 0.0 if y in (0.0, 1.0) else y * log(y) + (1 - y) * log(1 - y)


omega = 8000.0
tc = omega / R
t = 0.8 * tc


def g_order(eta):
    y = (1 - eta) / 2
    z = (1 + eta) / 2
    return -omega * (y * (1 - z) + (1 - y) * z) + R * t * (q(y) + q(z)) / 2


lo, hi = 1e-10, 1 - 1e-12
for _ in range(100):
    eta = (lo + hi) / 2
    derivative = -omega * eta + R * t * log((1 + eta) / (1 - eta)) / 2
    if derivative > 0:
        hi = eta
    else:
        lo = eta
eta = (lo + hi) / 2
assert g_order(eta) < g_order(0) and isclose(
    g_order(eta), g_order(-eta), abs_tol=1e-8
)
print(
    (
        f"Symmetric two-sublattice toy: Tc="
        f"{tc:.6f}"
        f" K; T/Tc=0.8, order parameter="
        f"{eta:.9f}"
    )
)
print(
    (
        f"Sites at x_B=0.5: y_B1="
        f"{(1 - eta) / 2:.9f}"
        f", y_B2="
        f"{(1 + eta) / 2:.9f}"
        f"; ordered-minus-disordered G="
        f"{g_order(eta) - g_order(0):.6f}"
        f" J/mol atoms"
    )
)
# Interstitial composition: vacancies count sites, not atoms.
yc = 0.1
nc = 3 * yc
nfe = 1.0
xc = nc / (nc + nfe)
sformula = -3 * R * q(yc)
print(
    (
        f"(Fe)1(C,Va)3 at y_C=0.1: x_C="
        f"{xc:.9f}"
        f"; entropy="
        f"{sformula:.6f}"
        f" J/(mol formula K)"
    )
)
print(
    (
        f"Same entropy per mole real atoms="
        f"{sformula / (1 + nc):.6f}"
        f"; vacancy-only sublattice contributes zero ideal entropy"
    )
)
# Ideal associated solution A+B <-> AB, total A=B=0.5 mol atoms.
keq = 9.0
z = (1 - sqrt(1 / (1 + keq))) / 2
na = nb = 0.5 - z
nab = z
nt = na + nb + nab
assert isclose((nab / nt) / ((na / nt) * (nb / nt)), keq, abs_tol=1e-12)
assert isclose(na + nab, 0.5) and isclose(nb + nab, 0.5)
print(
    (
        f"Association K=9: n_AB="
        f"{nab:.9f}"
        f", species fraction y_AB="
        f"{nab / nt:.9f}"
        f", overall atomic x_B=0.500000000"
    )
)
# Linearized observations for L(T)=A+B*(T-1000), with known independent
# Gaussian errors.
observations = [
    ("normalized enthalpy", [1.0, -1000.0], 9100.0, 100.0),
    ("L at 1000K", [1.0, 0.0], 5950.0, 200.0),
    ("L at 1500K", [1.0, 500.0], 4580.0, 200.0),
]
normal = [[0.0, 0.0], [0.0, 0.0]]
rhs = [0.0, 0.0]
for label, row, value, u in observations:
    for i in range(2):
        rhs[i] += row[i] * value / u**2
        for j in range(2):
            normal[i][j] += row[i] * row[j] / u**2
det = normal[0][0] * normal[1][1] - normal[0][1] ** 2
assert det > 0
cov = [
    [normal[1][1] / det, -normal[0][1] / det],
    [-normal[1][0] / det, normal[0][0] / det],
]
theta = [sum(cov[i][j] * rhs[j] for j in range(2)) for i in range(2)]
chi2 = sum(
    (sum(row[i] * theta[i] for i in range(2)) - value) ** 2 / u**2
    for _, row, value, u in observations
)
print(
    (
        f"Weighted fit: L(1000)="
        f"{theta[0]:.6f}"
        f" J/mol, slope="
        f"{theta[1]:.9f}"
        f" J/(mol K), chi2="
        f"{chi2:.9f}"
    )
)
print(
    (
        f"Known-noise covariance: u_L1000="
        f"{sqrt(cov[0][0]):.6f}"
        f", u_slope="
        f"{sqrt(cov[1][1]):.9f}"
        f", correlation="
        f"{cov[0][1] / sqrt(cov[0][0] * cov[1][1]):.9f}"
    )
)
for temp in (1000.0, 1500.0, 3000.0):
    row = [1.0, temp - 1000.0]
    pred = sum(row[i] * theta[i] for i in range(2))
    variance = sum(
        row[i] * cov[i][j] * row[j] for i in range(2) for j in range(2)
    )
    print(
        (
            f"Conditional model prediction T="
            f"{temp:.0f}"
            f" K: L="
            f"{pred:.6f}"
            f", standard uncertainty="
            f"{sqrt(variance):.6f}"
            f" J/mol"
        )
    )
# Check the exact weighted gradient, not only a printed optimizer status.
grad = [
    sum(
        row[i] * (sum(row[j] * theta[j] for j in range(2)) - value) / u**2
        for _, row, value, u in observations
    )
    for i in range(2)
]
assert max(abs(v) for v in grad) < 1e-10
# Every binary/ternary boundary misses an interior-only quaternary term.
hidden = -1e6 * (0.25**4)
assert -1e6 * 0 * 0.2 * 0.3 * 0.5 == 0.0
print(
    (
        f"Quaternary-only term: zero on every lower-order boundary, equimo"
        f"lar shift="
        f"{hidden:.6f}"
        f" J/mol atoms"
    )
)
print(
    (
        f"Five sublattices: 2 constituents each -> "
        f"{2 ** 5}"
        f" end members; 4 each -> "
        f"{4 ** 5}"
    )
)
print(
    (
        f"Four components: "
        f"{comb(4, 2)}"
        f" binaries, "
        f"{comb(4, 3)}"
        f" ternaries; subsystem counts do not prove extrapolation accuracy"
    )
)
known = {("FE", "VA"), ("CR", "VA"), ("FE", "C")}
all_members = {(a, b) for a in ("FE", "CR") for b in ("C", "VA")}
missing = sorted(all_members - known)
assert missing == [("CR", "C")]
print(
    (
        "Merged (Fe,Cr)1(C,Va)1: missing end member Cr:C must be explicitly a"
        "ssessed or labeled provisional"
    )
)
print(
    (
        "All model-constraint, weighted-gradient and balance checks PASS; no "
        "MCMC or commercial Ag-Cu reproduction claimed"
    )
)
