"""Original size, transport, electrochemical and database-scope audits.
Synthetic parameters only: not an assessed alloy/corrosion/process simulation.
"""
from math import pi, exp, log, log10, comb, isclose
R, F = 8.31446261815324, 96485.33212331002
vm, gamma, radius = 1e-5, 1., 10e-9
n = 4*pi*radius**3/(3*vm)
def area(moles):
    return 4*pi*(3*vm*moles/(4*pi))**(2/3)
average = gamma*area(n)/n
partial = 2*gamma*vm/radius
step = n*1e-5
finite = gamma*(area(n+step)-area(n-step))/(2*step)
assert isclose(average, 3000., rel_tol=1e-12)
assert isclose(partial, 2000., rel_tol=1e-12)
assert isclose(finite, partial, rel_tol=1e-9)
print(f'Spherical surface contribution, total/n and derivative: {average:.6f}, {partial:.6f} J/mol')
print(f'Centered finite-difference derivative: {finite:.6f} J/mol')
# Equal volumes and equal constrained radius, constant surface energies.
for rnm in (5., 10., 50.):
    shift = 3*vm*(.6-1.)/(rnm*1e-9)
    tm = 1000.+shift/10.
    assert abs(10*(1000-tm)+shift)<1e-9
    print(f'Fixed-radius two-state crossing, radius {rnm:g} nm: {tm:.6f} K')

# Toy parallel/series scalar paths; geometric weights are stipulated.
T, Q = 1000., 200000.
fgb, fp = .5e-9/(50e-9), 1e15*(.25e-9)**2
fb = 1-fgb-fp
assert min(fgb,fp,fb)>=0
mgb, mp = exp((1-.5)*Q/(R*T)), exp((1-.8)*Q/(R*T))
parallel = fb + fgb*mgb + fp*mp
series = 1/(fb + fgb/mgb + fp/mp)
assert series<=parallel and fb+fgb+fp==1.
print(f'Toy path fractions bulk/GB/pipe: {fb:.8f}, {fgb:.8f}, {fp:.8f}')
print(f'Toy effective/bulk mobility, parallel and series: {parallel:.6f}, {series:.6f}')

# Reduction Ox + mH+ + ne- = Red; all other activities unity.
T=298.15
slope=R*T*log(10)/F
print(f'Nernst decade factor at 298.15 K: {slope:.9f} V')
for ph in (0.,7.,14.):
    potential=.8-slope*ph  # m=n=2
    print(f'Illustrative m/n=1 reduction boundary, pH {ph:g}: {potential:.9f} V')
metal=-.2+slope/2*log10(1e-6)
print(f'Illustrative M2+/M boundary at ion activity 1e-6: {metal:.9f} V')

# M + .5 O2 -> MO, unit condensed activities and ideal oxygen gas.
T,dg0=1500.,-200000.
logp=2*dg0/(R*T*log(10))
assert abs(dg0-.5*R*T*log(10**logp))<1e-8
print(f'Toy oxide equilibrium log10(pO2/pstandard): {logp:.9f}')
print(f'Five-component binary/ternary subsystem counts: {comb(5,2)}, {comb(5,3)}')

# Two ternary models agree on every binary edge but differ in the interior.
for sign in (-1.,1.):
    interaction=sign*27000.*(1/3)**3
    assert sign*27000.*0.*.4*.6==0.
    print(f'Edge-invisible ternary interaction at center: {interaction:.6f} J/mol')
# Closed binary: pure A/B at zero, an optional AB line phase at x=.5.
z=.4
old=0.
sigma_fraction=2*min(z,1-z)
new=-500*sigma_fraction
assert new<=old
print(f'Missing-phase toy at x=.4, restricted/enlarged envelope: {old:.6f}, {new:.6f} J/mol')
# Stipulated bounded model discrepancy; not inferred confidence levels.
margin, error = -300., 500.
print(f'Energy-difference interval under stipulated bound: [{margin-error:.1f}, {margin+error:.1f}] J/mol')
