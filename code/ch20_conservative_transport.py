"""Original conservative transport benchmarks; no commercial solver or alloy fit."""
from math import exp,log,sqrt,pi,erf,erfc,sin,cos,ceil,isclose
R=8.31446261815324
length,diff,mean=2e-4,1e-13,.25
def exact_average(i,n,tau):
    a,b=i/n,(i+1)/n
    return mean+sum(-2/(m*pi)**2*exp(-(m*pi)**2*tau)*
        (sin(m*pi*b)-sin(m*pi*a))/(m*pi*(b-a))
        for m in range(1,301,2))
def tridiagonal(lower,diag,upper,rhs):
    n=len(diag);c=upper[:];d=rhs[:];pivot=diag[0]
    if n>1:c[0]/=pivot
    d[0]/=pivot
    for i in range(1,n):
        pivot=diag[i]-lower[i-1]*c[i-1]
        assert pivot>0
        if i<n-1:c[i]/=pivot
        d[i]=(d[i]-lower[i-1]*d[i-1])/pivot
    for i in range(n-2,-1,-1):d[i]-=c[i]*d[i+1]
    return d
def evolve(n,tau,variable=False):
    u=[.5*(i+.5)/n for i in range(n)]
    cell_d=[1. if not variable or i<n//2 else 4. for i in range(n)]
    face=[2*a*b/(a+b) for a,b in zip(cell_d,cell_d[1:])]
    steps=ceil(tau*n*n/.5);step=tau/steps;r=step*n*n
    lower=[-r*d for d in face];upper=lower[:]
    diag=[1+r*((face[i-1] if i else 0)+(face[i] if i<n-1 else 0)) for i in range(n)]
    mass0=sum(u)/n;var0=sum((x-mass0)**2 for x in u)/n
    max_drift=0.;largest_variance_rise=0.;oldvar=var0
    for _ in range(steps):
        u=tridiagonal(lower,diag,upper,u)
        mass=sum(u)/n;var=sum((x-mass0)**2 for x in u)/n
        max_drift=max(max_drift,abs(mass-mass0))
        largest_variance_rise=max(largest_variance_rise,var-oldvar)
        assert min(u)>-1e-12 and max(u)<.5+1e-12
        oldvar=var
    assert max_drift<1e-12 and largest_variance_rise<1e-14
    return u,steps,max_drift,oldvar
errors=[]
for n in (40,80,160):
    u,steps,drift,var=evolve(n,.1)
    err=max(abs(u[i]-exact_average(i,n,.1)) for i in range(n))
    errors.append(err)
    print(f'Closed diffusion N={n}: steps={steps}, max error={err:.9e}, mean drift={drift:.3e}, variance={var:.9e}')
assert errors[0]/errors[1]>3.8 and errors[1]/errors[2]>3.8
print(f'Convergence ratios={errors[0]/errors[1]:.6f},{errors[1]/errors[2]:.6f}; second-order combined refinement')
u,steps,drift,var=evolve(80,.1,True)
print(f'Variable-D closed specimen: mean={sum(u)/80:.12f}, range={min(u):.9f}..{max(u):.9f}; every-step conservation, positivity and variance checks PASS')
tol=.001
tau_est=log((2/pi**2)/tol)/pi**2
time_est=tau_est*length**2/diff
remainder=sum(2/(m*pi)**2*exp(-m*m*pi*pi*tau_est) for m in range(3,301,2))
assert remainder<1e-19
print(f'First-mode homogenization threshold {tol}: time={time_est/3600:.6f} h; higher-mode bound={remainder:.3e}')
# Spatially uniform prescribed temperature; no thermodiffusion term.
d0,q=1e-5,180000.
def diffusivity(t):return d0*exp(-q/(R*t))
def simpson(f,a,b,n):
    assert n%2==0
    h=(b-a)/n
    return h/3*(f(a)+f(b)+sum((4 if i%2 else 2)*f(a+i*h) for i in range(1,n)))
def temperature(t):return 1173+(1050-1173)*t/773 if t<=773 else 1050.
def clock(n):return simpson(lambda t:diffusivity(temperature(t)),0,773,n)+227*diffusivity(1050)
theta=clock(2000);coarse=clock(1000)
assert abs(theta-coarse)/theta<1e-9
print(f'Thermal diffusion clock={theta:.12e} m2; equivalent1050K time={theta/diffusivity(1050):.6f} s; quadrature relative difference={abs(theta-coarse)/theta:.3e}')
# One-sided planar local-equilibrium growth into a semi-infinite parent.
cs,ci,c0,di=.02,.4,.15,1e-12
ratio=(ci-c0)/(sqrt(pi)*(ci-cs))
lo,hi=0.,3.
for _ in range(100):
    lam=(lo+hi)/2
    if lam*exp(lam*lam)*erfc(lam)<ratio:lo=lam
    else:hi=lam
lam=(lo+hi)/2;t=1000.;position=2*lam*sqrt(di*t)
velocity=lam*sqrt(di/t)
flux=di*(ci-c0)*exp(-lam*lam)/(sqrt(pi*di*t)*erfc(lam))
assert isclose(velocity*(ci-cs),flux,rel_tol=1e-13)
excess=2*sqrt(di*t)*(ci-c0)*(exp(-lam*lam)/sqrt(pi)-lam*erfc(lam))/erfc(lam)
assert isclose(excess,(c0-cs)*position,rel_tol=1e-13)
print(f'One-sided interface: lambda={lam:.9f}, s1000={position*1e6:.6f} um, v={velocity:.12e} m/s; local jump and integrated solute PASS')
fraction=(ci-c0)/(ci-cs)
print(f'Finite closed equal-density equilibrium: solid fraction={fraction:.9f}; semi-infinite growth is not valid through finite saturation')
# Two coupled constant coefficients, stable single phase, independent mole fractions.
z=2*sqrt(1e-13*1000.)
eplus=erf(z/(2*sqrt(3e-13*1000.)));eminus=erf(z/(2*sqrt(1e-13*1000.)))
xa=.2+(.3*eplus+.1*eminus)/4
xb=.15+(.3*eplus-.1*eminus)/4
offset=xb-(.1+.5*(xa-.1))
assert 0<xa<1 and 0<xb<1 and xa+xb<1 and abs(offset)>.005
print(f'Coupled diffusion path at20um: xA={xa:.9f}, xB={xb:.9f}, deviation from end-member straight line={offset:.9f}')
# Additive interface and diffusion resistances, an explicitly defined toy law.
aa,ell=1e-14,1e-6
mixed=[]
for t in (1.,100.,10000.):
    s=2*aa*t/(sqrt(ell*ell+2*aa*t)+ell)
    exponent=(ell+s/2)/(ell+s)
    assert isclose(s*s/2+ell*s,aa*t,rel_tol=1e-13)
    mixed.append((t,s*1e6,exponent))
print('Mixed-control (time_s,position_um,local exponent): '+str([(t,round(s,9),round(p,9)) for t,s,p in mixed]))
rates=(10.,1.,.1,.001)
fourier=[4*1e-13*(100/r)/(50e-6*r**(-1/3))**2 for r in rates]
print('Back-diffusion scale (K/s, integratedD/halfspacing2): '+str([(r,round(f,9)) for r,f in zip(rates,fourier)]))
# Dilute continuous-growth partition law at prescribed fixed velocity.
ke,vd,initial,fs=.2,1.,.04,.9
partitions=[]
for speed in (0.,.1,1.,10.):
    k=(ke+speed/vd)/(1+speed/vd)
    liquid=initial*(1-fs)**(k-1)
    average_solid=(initial-(1-fs)*liquid)/fs
    assert isclose(fs*average_solid+(1-fs)*liquid,initial,abs_tol=1e-14)
    partitions.append((speed,k,liquid))
print('Fixed-speed trapping toy (v_m/s,k,remaining liquid x): '+str([(v,round(k,9),round(c,9)) for v,k,c in partitions]))
assert isclose((.9-ke)/(1-.9),7.)
print('All checks PASS; no commercial demo or multicomponent solidification history was reproduced')
