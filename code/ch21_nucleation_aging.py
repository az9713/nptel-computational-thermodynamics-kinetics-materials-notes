"""Synthetic nucleation, random growth, particle moments and strength checks."""
from math import pi,exp,expm1,log,log1p,sqrt,cos,isclose
KB=1.380649e-23;R=8.31446261815324
temp=623.15;gamma=.05;bulk=1e8;elastic=2e7;drive=bulk-elastic
radius=2*gamma/drive;barrier=16*pi*gamma**3/(3*drive**2)
def work(r):return 4*pi*gamma*r*r-4*pi*drive*r**3/3
assert isclose(work(radius),barrier,rel_tol=1e-14)
assert work(.9*radius)<barrier and work(1.1*radius)<barrier
shape=(2+cos(pi/3))*(1-cos(pi/3))**2/4
print(f'Capillarity: critical radius={radius*1e9:.9f} nm, zero crossing={3*gamma/drive*1e9:.9f} nm, barrier/kBT={barrier/(KB*temp):.9f}, flat60degree factor={shape:.9f}')
prefactor=1e32;rate=prefactor*exp(-barrier/(KB*temp));volume=1e-18;tobs=100.
prob=-expm1(-rate*volume*tobs)
print(f'Illustrative homogeneous events: rate={rate:.9e}/m3/s; probability in1e-18m3 over100s={prob:.9f}; mean wait={1/(rate*volume):.9f}s')
ii,uu=1e14,1e-8;k=pi*ii*uu**3/3
thresholds=[(-log1p(-x)/k)**.25 for x in (.01,.5,.99)]
print('Constant random nucleation/linear growth: K='+f'{k:.12e}/s4; t1%,t50%,t99%='+str([round(x,9) for x in thresholds]))
ns=1e18;k3=4*pi*ns*uu**3/3
print(f'Site-saturated linear growth: n=3, t50={(-log(.5)/k3)**(1/3):.9f}s; no added induction time')
# Exact moment update for ideal Poisson births and common linear radial growth,
# with I and u piecewise constant within a time interval.
def advance(m,i,u,dt):
    g=u*dt;a,b,c,d=m
    return (a+i*dt,b+g*a+i*u*dt**2/2,
        c+2*g*b+g*g*a+i*u*u*dt**3/3,
        d+3*g*c+3*g*g*b+g**3*a+i*u**3*dt**4/4)
base=(0.,0.,0.,0.)
ma=advance(advance(base,ii,uu,100),ii,uu,100)
mb=advance(advance(base,ii/8,uu*2,100),ii*8,uu/2,100)
xa=-expm1(-4*pi*ma[3]/3);xb=-expm1(-4*pi*mb[3]/3)
assert abs(xa-xb)>.05
assert isclose((ii/8)*(uu*2)**3,ii*uu**3)
assert isclose((ii*8)*(uu/2)**3,ii*uu**3)
print(f'Same isothermal K at both stages, different histories: X_A={xa:.9f}, X_B={xb:.9f}; TTT K alone does not determine CCT')
tm=1000.
def rates(t):
    assert 0<t<tm
    nuc=1e32*exp(-60000/(R*t)-1e9/(t*(tm-t)**2))
    dg=10000*(tm-t)/tm
    growth=1e-3*exp(-50000/(R*t))*(-expm1(-dg/(R*t)))
    return nuc,growth
def t_fraction(t,x):
    i,u=rates(t);return (-log1p(-x)/(pi*i*u**3/3))**.25
nose=min(range(300,951),key=lambda t:t_fraction(t,.01))
print(f'Synthetic TTT1% grid nose={nose}K, time={t_fraction(nose,.01):.9f}s; no real alloy calibration')
def cooling(cooling_rate,steps):
    hi,low=950.,300.;dt=(hi-low)/(cooling_rate*steps)
    m=(0.,0.,0.,0.);old_y=0.;found=None;target=-log1p(-.01)
    for j in range(steps):
        tt=hi-(j+.5)*(hi-low)/steps;i,u=rates(tt)
        m=advance(m,i,u,dt);y=4*pi*m[3]/3
        if found is None and y>=target:
            frac=(target-old_y)/(y-old_y)
            found=hi-(j+frac)*(hi-low)/steps
        old_y=y
    return found,-expm1(-old_y)
for speed in (1.,10.,100.):
    coarse=cooling(speed,2000);fine=cooling(speed,4000)
    assert coarse[0] is None and fine[0] is None or abs(coarse[0]-fine[0])<.01
    print(f'Synthetic continuous cooling {speed:g}K/s: T1%={fine[0]}, final X={fine[1]:.9f}; refinement temperature difference={None if fine[0] is None else abs(coarse[0]-fine[0])}')
# A two-bin geometric distribution; equal total molar density approximates volume=amount.
radii=(2e-9,8e-9);number=(2.5e21,2.5e21)
n=sum(number);rmean=sum(a*b for a,b in zip(number,radii))/n
f=4*pi/3*sum(a*b**3 for a,b in zip(number,radii))
wrong=4*pi*n*rmean**3/3
initial,cp,ceq=.0018,.25,.0001
matrix=(initial-f*cp)/(1-f);feq=(initial-ceq)/(cp-ceq)
assert 0<ceq<matrix<initial and 0<f<feq
print(f'Two-bin population: mean radius={rmean*1e9:.6f}nm, exact volume fraction={f:.12f}, mean-radius substitution={wrong:.12f}, ratio={f/wrong:.9f}')
print(f'Solute inventory: matrix x={matrix:.12f}; planar-equilibrium amount fraction={feq:.12f}; current fraction below inventory bound')
capgamma,vm=.1,1e-5;reservoir=1.12*ceq;d=1e-19
critical=2*capgamma*vm/(R*temp*log(reservoir/ceq))
growth=[]
for r in radii:
    surface=ceq*exp(2*capgamma*vm/(r*R*temp))
    dr=d/r*(reservoir-surface)/(cp-surface)
    growth.append(dr)
assert growth[0]<0<growth[1]
print(f'Gibbs-Thomson toy: reservoir critical radius={critical*1e9:.9f}nm; growth at2/8nm={growth[0]:.12e},{growth[1]:.12e}m/s')
combined=(sum(r**3 for r in radii))**(1/3)
area_ratio=combined**2/sum(r*r for r in radii)
assert area_ratio<1
print(f'Equal-volume area check: combined radius={combined*1e9:.9f}nm; new/old area={area_ratio:.9f}; geometric comparison, not coalescence simulation')
def strength(r):
    shear=100*sqrt(r/8e-9);bypass=100*8e-9/r
    return 50+3*min(shear,bypass)
print('Illustrative fixed-fraction strength (radius_nm,yield_MPa): '+str([(r,round(strength(r*1e-9),9)) for r in (2.,8.,32.)]))
print('All checks PASS; no KWN solver, commercial precipitation demo or experimental aging curve reproduced')
