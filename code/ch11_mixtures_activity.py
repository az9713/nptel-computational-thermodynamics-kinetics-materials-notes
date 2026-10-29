"""Declared bulk regular-solution and activity models; no fitted material data."""
from math import log, exp, lgamma, isclose
R=8.31446261815324
GA,GB,OMEGA=-10000.0,-20000.0,12000.0

def molar_g(x,t,omega=OMEGA):
    return (1-x)*GA+x*GB+R*t*((1-x)*log(1-x)+x*log(x))+omega*x*(1-x)

def chemical(x,t,omega=OMEGA):
    return GA+R*t*log(1-x)+omega*x*x,GB+R*t*log(x)+omega*(1-x)**2

def total_g(na,nb,t):
    n=na+nb
    return n*molar_g(nb/n,t)

def simpson(fn,a,b,n=1000):
    h=(b-a)/n
    return h*(fn(a)+fn(b)+sum((4 if j%2 else 2)*fn(a+j*h) for j in range(1,n)))/3

def main():
    t,x=800.0,0.3
    ma,mb=chemical(x,t)
    g=molar_g(x,t)
    assert isclose(g,(1-x)*ma+x*mb,abs_tol=1e-9)
    h=1e-5;na,nb=.7,.3
    da=(total_g(na+h,nb,t)-total_g(na-h,nb,t))/(2*h)
    db=(total_g(na,nb+h,t)-total_g(na,nb-h,t))/(2*h)
    assert isclose(da,ma,rel_tol=1e-9) and isclose(db,mb,rel_tol=1e-9)
    gd=(1-x)*(-R*t/(1-x)+2*OMEGA*x)+x*(R*t/x-2*OMEGA*(1-x))
    assert abs(gd)<1e-10
    print(f'800K xB0.3: g={g:.9f}; muA={ma:.9f}; muB={mb:.9f} J/mol; Euler and amount derivatives passed')
    aa=exp((ma-GA)/(R*t));ab=exp((mb-GB)/(R*t))
    gamma_a=aa/(1-x);gamma_b=ab/x
    print(f'Activities: aA={aa:.9f}; aB={ab:.9f}; gammaA={gamma_a:.9f}; gammaB={gamma_b:.9f}')
    pa,pb=100*aa,20*ab
    print(f'Matched-reference ideal vapor: pA={pa:.9f}; pB={pb:.9f} Pa')
    finite=total_g(.9,.3,t)-total_g(.7,.3,t)
    differential=ma*.2
    integral=simpson(lambda a:chemical(.3/(a+.3),t)[0],.7,.9)
    assert isclose(finite,integral,abs_tol=1e-8)
    print(f'Adding0.2molA: exact dG={finite:.9f}; initial-mu approximation={differential:.9f}; integrated mu={integral:.9f} J')
    tm=600.0;xmid=.5
    mix=molar_g(xmid,tm)-((1-xmid)*GA+xmid*GB)
    curvature=R*tm/(xmid*(1-xmid))-2*OMEGA
    # A symmetric split preserves mean composition and cancels the linear baseline.
    split=.5*(molar_g(.49,tm)+molar_g(.51,tm))-molar_g(.5,tm)
    assert mix<0 and curvature<0 and split<0
    print(f'600K midcomposition: mixing_g={mix:.9f} J/mol; curvature={curvature:.9f}; symmetric split minus homogeneous={split:.9f} J/mol')
    # Phase-specific standard potentials differ: equal physical mu does not require equal activity.
    ap=.4;standard_delta=3000.0;teq=1000.0
    aq=ap*exp(-standard_delta/(R*teq))
    assert isclose(R*teq*log(ap),standard_delta+R*teq*log(aq),abs_tol=1e-9)
    print(f'Equal-mu with3000J/mol standard offset at1000K: a_alpha={ap:.9f}; a_beta={aq:.9f}')
    sexact=lgamma(101)-lgamma(31)-lgamma(71)
    sasym=-100*(.3*log(.3)+.7*log(.7))
    print(f'100site30B count: exact S/k={sexact:.9f}; extensive Stirling term={sasym:.9f}')
    five=R*log(5)
    other=-R*(.4*log(.4)+4*.15*log(.15))
    print(f'Five-cation ideal entropy: equimolar={five:.9f}; fractions0.4,0.15x4={other:.9f} J/(mol_sites K); ideal T*S at1148.15K={1148.15*five:.9f} J/mol_sites')
    print('All Euler, Gibbs-Duhem, derivative, activity and finite-addition checks passed.')

if __name__=='__main__':main()
