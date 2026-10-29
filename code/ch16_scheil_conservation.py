"""Conservative synthetic Scheil limits; no assessed alloy or commercial GUI."""
from math import exp, log, isclose, sqrt

C0, K, CE = 0.01, 0.13, 0.12  # mass fractions, dimensionless partition ratio
TM, SLOPE = 933.0, -600.0  # kelvin; kelvin per unit solute mass fraction

def liquid(fs):
    return C0 * (1.0-fs)**(K-1.0)

def average_solid(fs):
    if fs == 0.0:
        return K*C0  # right-hand limit, no finite solid amount yet
    return C0*(1.0-(1.0-fs)**K)/fs

re = (CE/C0)**(1.0/(K-1.0)); fe = 1.0-re
te = TM+SLOPE*CE
tl = TM+SLOPE*C0
tsolid_eq = TM+SLOPE*C0/K
assert C0 < K*CE and 0 < re < 1
print(f'Synthetic path: liquidus={tl:.6f} K, equilibrium solidus={tsolid_eq:.6f} K, eutectic={te:.6f} K')
for fs in (0.5, 0.9):
    cl=liquid(fs); cs=K*cl; avg=average_solid(fs)
    ceq=C0/(K+(1-K)*(1-fs))
    assert isclose(fs*avg+(1-fs)*cl,C0,abs_tol=1e-14)
    print(f'fs={fs:.1f}: liquid={cl:.9f}, new solid={cs:.9f}, mean solid={avg:.9f}, equilibrium liquid={ceq:.9f}')
print(f'Eutectic starts: primary solid={fe:.9f}, residual liquid={re:.9f}, liquid composition={liquid(fe):.9f}')
cbeta=0.98; calpha=K*CE
beta_in_eut=(CE-calpha)/(cbeta-calpha)
final_beta=re*beta_in_eut
primary_solute=fe*average_solid(fe)
assert isclose(primary_solute+re*CE,C0,abs_tol=1e-14)
assert isclose(primary_solute+re*(1-beta_in_eut)*calpha+final_beta*cbeta,C0,abs_tol=1e-14)
print(f'Terminal mixture: beta within eutectic={beta_in_eut:.9f}, final beta mass fraction={final_beta:.9f}, final total solute={C0:.9f}')

def stepped(n):
    r=1.0; old=C0; stored=0.0
    for j in range(1,n+1):
        new=C0+(CE-C0)*j/n; cs=K*new
        theta=(new-old)/(new-cs)
        assert 0 <= theta <= 1
        df=r*theta; rnew=r*(1-theta)
        stored+=df*cs
        assert isclose(stored+rnew*new,C0,abs_tol=2e-14)
        r,old=rnew,new
    return r,stored
errors=[]
for n in (20,40,80,160):
    rn,mn=stepped(n); err=abs(rn-re);errors.append(err)
    print(f'Liquid-composition steps={n}: residual={rn:.9f}, absolute error={err:.9f}, conservation PASS')
assert all(b<a for a,b in zip(errors,errors[1:]))
# Composition nonlinear enthalpy is not preserved by replacing all layers by their mean.
first=C0*(1-re**K)
second=(K*C0)**2*(1-re**(2*K-1))/(2*K-1)
variance=second/fe-(first/fe)**2
assert variance > 0
B=1e6  # J/kg coefficient multiplying c^2, purely illustrative
print(f'Frozen-solid concentration variance={variance:.12f}; layer-minus-mean enthalpy={B*variance:.6f} J/kg of primary solid')
# Different time histories give exactly the same classical path.
ls=1e-5; ds=1e-14; dm=1e-8; mix_length=1e-5
for time in (0.1,100.0):
    fos=ds*time/ls**2; fol=dm*time/mix_length**2
    print(f'Clock illustration t={time:.1f} s: Fo_s={fos:.6f}, Fo_liquid={fol:.1f}; path itself has no clock')
assert isclose(liquid(0),C0) and isclose(average_solid(0),K*C0)
for k in (1.0,2.0):
    cl=C0*(1-0.5)**(k-1)
    assert (isclose(cl,C0) if k==1 else cl<C0)
print('Endpoint/partition checks: k=1 no segregation; k>1 depletes liquid; classical path stopped before unphysical divergence')
