"""Independent quadrature and caloric-model checks; no fitted material data."""
from math import exp, expm1, sqrt, pi, isclose
KB=1.380649e-23
NA=6.02214076e23
R=KB*NA

def simpson(f,a,b,n=4000):
    assert n > 0 and n%2 == 0
    step=(b-a)/n
    return step/3*(f(a)+f(b)+4*sum(f(a+i*step) for i in range(1,n,2))+2*sum(f(a+i*step) for i in range(2,n,2)))

def scaled_speed_density(s):
    return 4/sqrt(pi)*s*s*exp(-s*s)

def einstein_ratio(reduced_temperature):
    x=1/reduced_temperature
    # Equivalent to x² exp(x)/(exp(x)-1)² without positive exponential overflow.
    return x*x*exp(-x)/(-expm1(-x))**2

def cp(temperature):
    return 20.0+0.01*temperature+10000.0/temperature**2

def delta_h(t1,t2):
    return 20*(t2-t1)+0.005*(t2*t2-t1*t1)+10000*(1/t1-1/t2)

def main():
    norm=simpson(scaled_speed_density,0,10)
    first=simpson(lambda s:s*scaled_speed_density(s),0,10)
    second=simpson(lambda s:s*s*scaled_speed_density(s),0,10)
    assert isclose(norm,1,abs_tol=1e-10)
    assert isclose(first,2/sqrt(pi),abs_tol=1e-10)
    assert isclose(second,1.5,abs_tol=1e-10)
    temperature,molar_mass=300.,0.020 # stipulated 20 g/mol, not tabulated natural neon
    particle_mass=molar_mass/NA
    v0=sqrt(2*KB*temperature/particle_mass)
    speeds=[v0,first*v0,sqrt(second)*v0]
    assert speeds[0]<speeds[1]<speeds[2]
    volume=R*temperature/100000
    total_u=1.5*R*temperature
    assert isclose(1.5*100000*volume,total_u)
    # Same T and p, twice the amount: extensive variables double; molar volume does not.
    assert isclose((2*R*temperature/100000)/2,volume)
    analytic=delta_h(300,900)
    assert isclose(analytic,simpson(cp,300,900),rel_tol=1e-11)
    assert isclose(einstein_ratio(100),1,rel_tol=1e-4)
    for reduced in [.1,.5,1,2,10]:assert 0<einstein_ratio(reduced)<1
    eps=1e-3
    assert isclose((delta_h(300,600+eps)-delta_h(300,600-eps))/(2*eps),cp(600),rel_tol=1e-9)
    print('Exact SI product R = %.12f J/(mol K).' % R)
    print('Scaled Maxwell integrals: norm %.9f; mean %.9f; second %.9f.' % (norm,first,second))
    print('At 300 K, stipulated 20 g/mol: most-probable / mean / rms = ' + ' / '.join('%.6f' % v for v in speeds)+' m/s.')
    print('One mole at 300 K, 100 kPa: V = %.9f m3; monatomic U = %.6f J.' % (volume,total_u))
    print('Monatomic Cv,m = %.9f; Cp,m = %.9f; difference = %.9f J/(mol K).' % (1.5*R,2.5*R,R))
    print('Synthetic Cp integral, 300 -> 900 K: Delta H,m = %.9f J/mol.' % analytic)
    print('Einstein Cv/(3R) at T/Theta = 0.1, 0.5, 1, 2, 10: '+', '.join('%.9f' % einstein_ratio(t) for t in [.1,.5,1,2,10]))
    print('Checks: independent speed moments, scaling, caloric integral/derivative and quantum limits PASS.')

if __name__=='__main__':main()
