"""Synthetic ideal solid/liquid equilibrium and exact eutectic balances."""
from math import exp, expm1, log, isclose
from fractions import Fraction as F

R, TMA, TMB, HA, HB = 8.314462618, 1000.0, 1500.0, 10000.0, 12000.0


def endpoints(T):
    # Delta g_i = H_i*(1-T/T_mi), K_i = x_i^s/x_i^l.
    a, b = HA/R*(1/T-1/TMA), HB/R*(1/T-1/TMB)
    ka, kb = exp(a), exp(b)
    xl = -expm1(a)/(kb-ka)
    return xl, kb*xl


def invert_endpoint(x0, which):
    lo, hi = TMA, TMB
    for _ in range(70):
        mid = (lo+hi)/2
        if endpoints(mid)[which] < x0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def mix(x):
    return sum(z*log(z) for z in (x, 1-x) if z)


def free_energy(x, T, liquid=False):
    base = ((1-x)*HA*(1-T/TMA)+x*HB*(1-T/TMB)) if liquid else 0
    return base+R*T*mix(x)


if __name__ == '__main__':
    T, x0 = 1200.0, 0.45
    xl, xs = endpoints(T)
    assert 0 < xl < x0 < xs < 1
    fs = (x0-xl)/(xs-xl)
    assert isclose((1-fs)*xl+fs*xs, x0, abs_tol=1e-14)
    # Equal tangent slope and intercept in J/mol and J/mol per atomic fraction.
    slope = R*T*log(xs/(1-xs))
    liquid_slope = HB*(1-T/TMB)-HA*(1-T/TMA)+R*T*log(xl/(1-xl))
    intercept = free_energy(xs,T)-xs*slope
    liquid_intercept = free_energy(xl,T,True)-xl*liquid_slope
    assert abs(slope-liquid_slope) < 1e-9
    assert abs(intercept-liquid_intercept) < 1e-9
    minimum_gap = min(free_energy(i/5000,T,p)-(intercept+i/5000*slope)
                      for i in range(5001) for p in (False,True))
    assert minimum_gap >= -1e-9
    for i in range(101):
        z=xl+(xs-xl)*i/100
        f=(z-xl)/(xs-xl)
        assert -1e-12 <= f <= 1+1e-12
        assert abs((1-f)*xl+f*xs-z) < 1e-12
    ts, tl = invert_endpoint(x0,1), invert_endpoint(x0,0)
    assert ts < T < tl
    print(f'At T={T:.1f} K: x_B(liquid)={xl:.9f}, x_B(solid)={xs:.9f}')
    print(f'At overall x_B={x0:.2f}: solid atomic-amount fraction={fs:.9f}')
    print(f'Solidus={ts:.6f} K; liquidus={tl:.6f} K')
    print('101 lever balances and 10002 sampled support-line inequalities passed.')
    print(f'Minimum sampled tangent gap={minimum_gap:.9g} J/mol')
    # Stipulated eutectic mass-fraction endpoints, independent of ideal model.
    ca, ce, cb, overall = F(1,10), F(2,5), F(9,10), F(1,4)
    liquid = (overall-ca)/(ce-ca)
    eut_beta = (ce-ca)/(cb-ca)
    beta_total = liquid*eut_beta
    assert liquid == F(1,2) and eut_beta == F(3,8)
    assert beta_total == (overall-ca)/(cb-ca) == F(3,16)
    print('Hypoeutectic construction: remaining liquid=',liquid,
          '; beta inside eutectic=',eut_beta,'; total beta=',beta_total)
    print('All numbers are synthetic teaching models, not assessed alloy data.')