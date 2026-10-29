"""Original path/reference/T-zero audits with executed pycalphad0.11.1."""
from pathlib import Path
from math import exp,log,isclose
import numpy as np
import pycalphad
from pycalphad import Database,calculate,equilibrium,variables as v
assert pycalphad.__version__=='0.11.1'
R=float(v.R);db=Database(str(Path('data/ch15_tzero.tdb')))
# Independent linear conditions in Co,Cr,Fe,Ni,V order.
A=np.array([[1,1,1,1,1],[0,0,0,0,1],[1,-2,0,0,0],
            [0,1,-1,0,0],[0,0,1,-1,0]],dtype=float)
x=np.linalg.solve(A,np.array([1,.2,0,0,0]))
assert np.linalg.matrix_rank(A)==5 and np.allclose(x,[.32,.16,.16,.16,.2])
print('Five-component path: rank 5, fractions (.32,.16,.16,.16,.20) at xV=.20')

# Tangent geometry requires matching the ordinate and composition measure.
MA,MB=50.,100.;slope=2000.;intercept=-1000.
for xx in [.2,.8]:
    w=xx*MB/((1-xx)*MA+xx*MB)
    dxdw=2/(2-w)**2
    gm=slope*xx+intercept
    gmass=gm/((1-xx)*MA+xx*MB)
    assert isclose(gmass,-20+30*w,abs_tol=1e-12)
    print(f'Tangent contact x={xx:.1f}: w={w:.9f}, molar-plot slope={slope*dxdw:.6f}, specific g={gmass:.9f} J/g')
print('Correct mass-specific supporting line: -20 + 30w J/g')

def energy(phase,T,xx,output='GM'):
    res=calculate(db,['CU','NI'],[phase],T=T,P=100000,N=1,
                  points=np.array([[1-xx,xx]]),output=output)
    return float(res[output].values.squeeze())

z=.4;T0=1400-800*z
p=energy('PSOL',T0,z);q=energy('QSOL',T0,z)
assert abs(p-q)<1e-6
print(f'T0 at xNi=.4: {T0:.6f} K; equal homogeneous GM={p:.4f} J/mol')
D=4000.;delta=10*(1000-T0)
u=exp((D+delta)/(R*T0));vv=exp((delta-D)/(R*T0))
b=(u-1)/(u-vv);a=vv*b
muA=R*T0*log(1-a);m=D+R*T0*log(a/(1-a))
f=(z-a)/(b-a)
res=equilibrium(db,['CU','NI'],['PSOL','QSOL'],
                {v.T:T0,v.P:100000,v.N:1,v.X('NI'):z},calc_opts={'pdens':1000})
records=sorted((str(n),float(xx),float(ff)) for n,xx,ff in zip(
    res.Phase.values.reshape(-1),res.X.sel(component='NI').values.reshape(-1),res.NP.values.reshape(-1))
    if n and np.isfinite(ff) and ff>1e-9)
assert [r[0] for r in records]==['PSOL','QSOL']
assert abs(records[0][1]-a)<1e-7 and abs(records[1][1]-b)<1e-7
assert abs(records[1][2]-f)<1e-7
assert abs(sum(xx*ff for _,xx,ff in records)-z)<1e-9
assert abs(float(res.GM.values.squeeze())-(muA+m*z))<1e-6
print(f'Coexistence at T0: xP={a:.6f}, xQ={b:.6f}, fQ={f:.6f}, GM={float(res.GM.values.squeeze()):.4f}')
print(f'Homogeneous equality lies {p-float(res.GM.values.squeeze()):.4f} J/mol above global equilibrium')
# Positive curvature makes these common tangents global for the two branches.
assert a<z<b and p>float(res.GM.values.squeeze())

T=1000.
raw=energy('QSOL',T,z)-energy('PSOL',T,z)
wrong=energy('QSOL',T,z,'GM_MIX')-energy('PSOL',T,z,'GM_MIX')
assert abs(raw-800)<1e-6 and abs(wrong)<1e-6
wrong=0.0 if abs(wrong)<1e-6 else wrong
print(f'Common reference raw gQ-gP={raw:.4f}; independent phase mixing references difference={wrong:.4f}')
# A common P-endmember subtraction preserves the phase difference.
assert isclose((energy('QSOL',T,z)-4000*z)-(energy('PSOL',T,z)-4000*z),raw,abs_tol=1e-8)
print('Common elemental reference subtraction preserves phase competition: PASS')
# Parent Q to product P at the same z; below T0, the bulk chemical change favors P.
for T in [1000.,1080.,1160.]:
    dg=energy('PSOL',T,z)-energy('QSOL',T,z)
    assert abs(dg-10*(T-T0))<1e-6
    dg=0.0 if abs(dg)<1e-6 else dg
    print(f'Partitionless Q to P, T={T:.0f}: delta g={dg:.4f} J/mol')

# Root localization needs more than one sign test over a wide interval.
rootfn=lambda T:(T-800)*(T-1200)/1000
assert rootfn(600)>0 and rootfn(1400)>0
assert rootfn(1000)<0 and rootfn(800)==rootfn(1200)==0
print('Two-root counterexample: endpoints600/1400 positive; roots800 and1200; interior1000 negative')
touch=lambda T:(T-1000)**2/1000
assert touch(900)>0 and touch(1100)>0 and touch(1000)==0
print('Tangential root: 1000 K, no sign change and no stability-order reversal')
