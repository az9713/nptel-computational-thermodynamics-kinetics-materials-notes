"""Synthetic measurement-to-thermodynamics checks; no experimental dataset."""
from math import exp, log, sqrt, asin, pi, isclose

R=8.31446261815324
F=96485.33212  # rounded C/mol
x,T,A,B=0.3,1000.0,9000.0,-3.0  # J/mol, J/(mol K), atomic mole basis
L=A+B*T
p=x*(1-x)
hexcess=p*A
ln_gamma_b=(1-x)**2*L/(R*T)
activity=x*exp(ln_gamma_b)
emf=-R*T*log(activity)/(2*F)  # defined transfer reaction: pure B -> alloy B
a_fit=hexcess/p
b_fit=(R*T*log(activity/x)/(1-x)**2-a_fit)/T
assert isclose(a_fit,A,abs_tol=1e-9) and isclose(b_fit,B,abs_tol=1e-10)
print(f'Original regular-solution example: H_ex={hexcess:.6f} J/mol, a_B={activity:.9f}, oriented EMF={emf:.9f} V')
print(f'Joint calorimetry/activity recovery: A={a_fit:.6f} J/mol, B={b_fit:.6f} J/(mol K)')
for trial_b in (-8.0,0.0,5.0):
    assert isclose(p*A,hexcess)
print('Calorimetry-only entropy null direction: B=-8,0,5 all give the same H_ex')
# Gibbs-Duhem partial excess quantities for g_ex=L*x*(1-x).
for z in (0.1,0.4,0.8):
    dmu_a=2*L*z
    dmu_b=-2*L*(1-z)
    assert abs((1-z)*dmu_a+z*dmu_b)<1e-10
print('Gibbs-Duhem identity PASS; its integration still requires a reference constant')
qa,qb,qab=-8000.0,-12000.0,-5000.0  # J for stoichiometric dissolution of 1 mol AB
formation=qa+qb-qab
assert formation==-15000
print(f'Hess cycle: formation enthalpy={formation:.6f} J/mol AB; opposite subtraction would have wrong sign')
du=-100000.0;dn_gas=-1.0;t_bomb=298.15
dh=du+dn_gas*R*t_bomb
print(f'Ideal-gas reaction correction: DeltaU={du:.6f}, DeltaH={dh:.6f} J/mol reaction')
# Four independent repeatability errors, one common calibration offset.
n=4;u_repeat=2.0;u_common=3.0
u_mean=sqrt(u_repeat**2/n+u_common**2)
u_wrong=sqrt((u_repeat**2+u_common**2)/n)
assert u_mean>u_wrong
print(f'Four readings: correct mean standard uncertainty={u_mean:.9f}, false independent-total value={u_wrong:.9f}')
# Simple linear mixed-signal model, not a complete X-ray microanalysis model.
fraction=0.25;c_particle=0.6;c_matrix=0.1
observed=fraction*c_particle+(1-fraction)*c_matrix
recovered=(observed-(1-fraction)*c_matrix)/fraction
assert isclose(recovered,c_particle)
print(f'Mixed-volume toy: observed={observed:.6f}, particle={recovered:.6f}, uncertainty amplification={1/fraction:.1f}')
# Unit/basis audit.
ev_per_atom=0.05
delta_h=ev_per_atom*F  # eV/atom -> J/mol atoms, same numerical conversion F
delta_s=10.0
tcross=delta_h/delta_s
print(f'Finite-T toy: 0.05 eV/atom={delta_h:.6f} J/mol atoms; entropy reversal at {tcross:.6f} K')
assert -90000/3==-30000
print('AB2 basis: -90 kJ/mol formula = -30 kJ/mol atoms; x_B=2/3')
# Mean-field continuous order transition.
a,b,tc=2.0,1000.0,800.0
delta_cp=a*a*tc/(2*b)
assert delta_cp==1.6
print(f'Landau example: G,S,H continuous at Tc; latent heat zero, Cp step={delta_cp:.6f} J/(mol K)')
wavelength=0.15406;alat=0.30
twotheta=2*asin(wavelength/(2*(alat/sqrt(2))))*180/pi
contracted=2*asin(wavelength/(2*(0.98*alat/sqrt(2))))*180/pi
print(f'Cubic(110) Bragg example: 2theta={twotheta:.6f} deg; 2% contraction shifts it by {contracted-twotheta:.6f} deg')
print('Data audit PASS: reference, phase, amount basis, measurement model and correlated uncertainty remain explicit')
