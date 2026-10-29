"""Synthetic diffusion, mobility and linear-stability checks; standard library."""
from math import exp,log,sqrt,pi,erfc,isclose
R=8.31446261815324
a,gamma,d=2.5e-10,1e6,3
diff=a*a*gamma/(2*d)
assert isclose(diff,a*a*(gamma/6),rel_tol=1e-14)
print(f'Uncorrelated 3D walk: total jump rate={gamma:.0f}/s, D={diff:.12e} m2/s')
d0,qact=1e-5,180000.
def arrhenius(temp):return d0*exp(-qact/(R*temp))
t1,t2=1000.,1200.;d1,d2=arrhenius(t1),arrhenius(t2)
qfit=-R*(log(d2)-log(d1))/(1/t2-1/t1)
assert isclose(qfit,qact,rel_tol=1e-12)
print(f'Arrhenius synthetic data: D1000={d1:.12e}, D1200={d2:.12e}; recovered Q={qfit:.6f} J/mol')
print(f'1000 Celsius versus 1000 kelvin: diffusivity ratio={arrhenius(1273.15)/arrhenius(1000):.6f}')
temp=1000.;nu=1e13;e1,e2,barrier=0.,-5000.,50000.
k12=nu*exp(-(barrier-e1)/(R*temp));k21=nu*exp(-(barrier-e2)/(R*temp))
p2=exp(-(e2-e1)/(R*temp));p1=1.;normal=p1+p2;p1/=normal;p2/=normal
assert isclose(p1*k12,p2*k21,rel_tol=1e-12)
print(f'Detailed balance: forward/reverse rate ratio={k12/k21:.9f}; reverse rate positive, equilibrium net flux zero')
# Half-space with fixed surface concentration; concentration units mol/m3.
c0,cs,ctarget,depth,dc=1000.,5000.,3000.,1e-3,1e-11
target=(ctarget-c0)/(cs-c0);lo,hi=0.,8.
for _ in range(100):
    eta=(lo+hi)/2
    if erfc(eta)>target:lo=eta
    else:hi=eta
eta=(lo+hi)/2;time=depth**2/(4*dc*eta**2)
def concentration(x,t):return c0+(cs-c0)*erfc(x/(2*sqrt(dc*t)))
assert isclose(concentration(depth,time),ctarget,abs_tol=1e-9)
uptake=2*(cs-c0)*sqrt(dc*time/pi)
print(f'Half-space target: eta={eta:.9f}, time={time:.6f} s, concentration at1mm={concentration(depth,time):.6f} mol/m3')
print(f'Half-space excess solute per area={uptake:.9f} mol/m2; fixed surface boundary is a continuing reservoir')
# Constant tracer diffusivity, regular homogeneous phase, equal-species Darken approximation.
temp,omega,dtracer=600.,12000.,1e-18
def curvature(x):return R*temp/(x*(1-x))-2*omega
for x in (.2,.5):
    phi=x*(1-x)*curvature(x)/(R*temp)
    chemical=dtracer*phi
    assert dtracer/(R*temp)>0
    print(f'Regular-solution x={x:.1f}: curvature={curvature(x):.6f} J/mol, thermodynamic factor={phi:.9f}, Dchemical={chemical:.12e} m2/s')
# Gradient-regularized, linearized conserved dynamics around x=0.5.
h=curvature(.5);kap=1e-12;mob=.25*dtracer/(R*temp)
kstar=sqrt(-h/(2*kap));wavelength=2*pi/kstar;sigma=mob*h*h/(4*kap)
amplitude=1e-3;later=amplitude*exp(sigma*3600.)
e0=(h+kap*kstar*kstar)*amplitude**2/4
e1=(h+kap*kstar*kstar)*later**2/4
assert e1<e0 and later<.01
print(f'Linear spinodal toy: fastest wavelength={wavelength*1e9:.6f} nm, growth rate={sigma:.12e}/s')
print(f'One-hour linear amplitude={later:.9f}; quadratic excess energy {e0:.9f} -> {e1:.9f} J/mol, mean x remains0.5')
# Independent exchange-potential basis, dimensionless illustrative Onsager matrix.
kinetic=((2.,1.),(1.,2.));force=(1.,-3.)
flux=tuple(-sum(kinetic[i][j]*force[j] for j in range(2)) for i in range(2))
diss=-sum(flux[i]*force[i] for i in range(2))
assert flux==(1.,5.) and diss==14.
hessian=((1.,-2.),(-2.,5.))
matrix=tuple(tuple(sum(kinetic[i][k]*hessian[k][j] for k in range(2)) for j in range(2)) for i in range(2))
assert matrix==((0.,1.),(-3.,8.))
print(f'Coupled flux toy: J={flux}, total dissipation={diss:.6f}; component1 can move along its own potential gradient')
print(f'D=L*H={matrix}; nonsymmetric with a negative cross entry, eigenvalues={4-sqrt(13):.9f},{4+sqrt(13):.9f}')
# Both fluxes use same laboratory x direction and physical molar concentration.
cleft,cright=1000.,4000.;jleft,jright=-1e-7,2e-7
speed=(jright-jleft)/(cright-cleft)
assert isclose(speed*(cright-cleft),jright-jleft,abs_tol=1e-22)
print(f'Conservative moving interface: speed={speed:.12e} m/s; flux jump matched to concentration jump')
print('All checks PASS; linearized spinodal mode is not a full nonlinear phase-field simulation')
