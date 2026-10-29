"""Exact teaching constructions: composition, sublattices and phase balances."""
from fractions import Fraction as F
from itertools import combinations
from math import comb

# Eight labelled sites: first four alpha, last four beta; exactly four A.
hist = {k: 0 for k in range(5)}
for occupied in combinations(range(8), 4):
    hist[sum(i < 4 for i in occupied)] += 1
assert list(hist.values()) == [comb(4, k)**2 for k in range(5)]
assert sum(hist.values()) == comb(8, 4) == 70
print('A on alpha | eta | labelled configurations')
for k, count in hist.items():
    print(k, F(k, 2)-1, count, sep=' | ')
# All configurations have x_A=1/2, but eta ranges from -1 to 1.
assert all(F(k + (4-k), 8) == F(1, 2) for k in hist)

# Site fractions count vacancies; atomic fractions do not.
x_B, y_B = F(20, 80), F(20, 100)
assert x_B == F(1, 4) and y_B == F(1, 5)
print('Vacancy example: atomic x_B =', x_B, '; site y_B =', y_B)
# Generic interstitial model: r=1 interstitial site per host atom.
r, y = F(1), F(1, 4)
x = r*y/(1+r*y)
assert x == F(1, 5) and x/(r*(1-x)) == y
print('Interstitial example: y_B =', y, '; atomic x_B =', x)

# Explicitly rounded molar masses, g/mol: Fe 56, C 12.
x_C, M_C, M_Fe = F(1, 4), F(12), F(56)
w_C = x_C*M_C/(x_C*M_C + (1-x_C)*M_Fe)
assert w_C == F(1, 15)
print('Fe3C rounded-mass construction: carbon mass fraction =', w_C)

# Original consistent two-phase table: alpha (18A,26B), eta (12A,4B).
counts = [(18, 26), (12, 4)]
N = sum(a+b for a, b in counts)
phase_f = [F(a+b, N) for a, b in counts]
phase_x = [F(a, a+b) for a, b in counts]
x_A = sum(f*c for f, c in zip(phase_f, phase_x))
assert phase_f == [F(11, 15), F(4, 15)]
assert phase_x == [F(9, 22), F(3, 4)] and x_A == F(1, 2)
lever_gamma = (x_A-phase_x[0])/(phase_x[1]-phase_x[0])
assert lever_gamma == phase_f[1]
print('Phase atom-amount fractions:', *phase_f)
print('Phase atomic compositions x_A:', *phase_x)
print('Overall x_A =', x_A, '; lever gamma fraction =', lever_gamma)
print('All exact conservation and enumeration checks passed.')