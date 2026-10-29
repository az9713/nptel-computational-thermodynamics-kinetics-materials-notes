"""Declared caloric and sensor models, not a real-alloy assessment."""
from math import exp, log, isclose
CP, LATENT, MASS, POWER = 500.0, 100000.0, 1.0, 100.0

def pure(t):
    h = 250000.0 - POWER*t/MASS
    if h > LATENT:
        return 1000.0 + (h-LATENT)/CP, 0.0, h
    if h >= 0:
        return 1000.0, 1-h/LATENT, h
    return 1000.0 + h/CP, 1.0, h

def mushy(t):
    h = 300000.0 - POWER*t/MASS
    if h > 200000:
        return 900.0 + (h-LATENT)/CP, 0.0, h
    if h >= 0:
        temperature = 900.0 + h/1000.0
        return temperature, (1100-temperature)/200, h
    return 900.0 + h/CP, 1.0, h

def sensor_step(start, slope, duration, observed, tau):
    """Exact response to a linear true-temperature segment."""
    end = start + slope*duration
    return end-slope*tau+(observed-start+slope*tau)*exp(-duration/tau)

def relaxation(t, equilibrium=0.8, tau=100.0):
    return equilibrium*(1-exp(-t/tau))

def main():
    for t in [0,1500,2000,2500,3500]:
        temperature, solid, h = pure(t)
        assert isclose(h,CP*(temperature-1000)+LATENT*(1-solid),abs_tol=1e-8)
    for t in range(0,4001):
        temperature, solid, h = mushy(t)
        assert isclose(h,CP*(temperature-900)+LATENT*(1-solid),abs_tol=1e-8)
        assert 0 <= solid <= 1
    assert pure(1500)[0] == pure(2500)[0] == 1000
    assert isclose(pure(2000)[1],0.5)
    assert mushy(1000)[0] == 1100 and mushy(3000)[0] == 900
    assert isclose((mushy(2001)[0]-mushy(2000)[0]),-0.1,abs_tol=1e-10)
    # On a ramp cooling at 0.2 K/s, the initially equilibrated sensor reads high.
    observed=sensor_step(1300,-0.2,100,1300,10)
    offset=observed-1280
    assert isclose(offset,2*(1-exp(-10)),abs_tol=1e-10)
    # Equal final enthalpy in an approximately adiabatic recalescence burst.
    undercool=20.0
    burst_fraction=CP*undercool/LATENT
    assert isclose(LATENT-CP*undercool,LATENT*(1-burst_fraction))
    assert isclose(relaxation(100*log(20))/0.8,0.95)
    print('Pure: liquidus arrival 1500 s; arrest 1000 s; completion 2500 s.')
    print('Pure at 2000 s: T = %.1f K; solid mass fraction = %.3f.' % pure(2000)[:2])
    print('Range model: first solid 1000 s; last liquid 3000 s; slope -0.100 K/s.')
    print('Sensor at 100 s: true 1280.000 K; observed %.6f K; offset +%.6f K.' % (observed,offset))
    print('Recalescence: 20 K undercooling -> %.3f solid mass fraction.' % burst_fraction)
    print('Relaxation at 10, 100, 300 s: ' + ', '.join('%.6f' % relaxation(t) for t in [10,100,300]))
    print('95%% of equilibrium amount: %.6f s.' % (100*log(20)))
    print('Diffusion estimates: %.0f s (%.6f days); %.0f s.' % (1e-8/1e-14,1e-8/1e-14/86400,1e-8/1e-9))
    print('Checks: 4001 caloric identities, pure-phase boundaries, sensor sign and burst enthalpy PASS.')

if __name__ == '__main__':
    main()
