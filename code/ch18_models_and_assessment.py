"""Original CEF, association and weighted-assessment examples; no assessed alloy."""
from math import log,sqrt,isclose,comb
R=8.31446261815324

def q(y):
    return 0.0 if y in (0.0,1.0) else y*log(y)+(1-y)*log(1-y)

omega=8000.0;tc=omega/R;t=.8*tc
def g_order(eta):
    y=(1-eta)/2;z=(1+eta)/2
    return -omega*(y*(1-z)+(1-y)*z)+R*t*(q(y)+q(z))/2
lo,hi=1e-10,1-1e-12
for _ in range(100):
    eta=(lo+hi)/2
    derivative=-omega*eta+R*t*log((1+eta)/(1-eta))/2
    if derivative>0:hi=eta
    else:lo=eta
eta=(lo+hi)/2
assert g_order(eta)<g_order(0) and isclose(g_order(eta),g_order(-eta),abs_tol=1e-8)
print(f'Symmetric two-sublattice toy: Tc={tc:.6f} K; T/Tc=0.8, order parameter={eta:.9f}')
print(f'Sites at x_B=0.5: y_B1={(1-eta)/2:.9f}, y_B2={(1+eta)/2:.9f}; ordered-minus-disordered G={g_order(eta)-g_order(0):.6f} J/mol atoms')
# Interstitial composition: vacancies count sites, not atoms.
yc=.1;nc=3*yc;nfe=1.
xc=nc/(nc+nfe);sformula=-3*R*q(yc)
print(f'(Fe)1(C,Va)3 at y_C=0.1: x_C={xc:.9f}; entropy={sformula:.6f} J/(mol formula K)')
print(f'Same entropy per mole real atoms={sformula/(1+nc):.6f}; vacancy-only sublattice contributes zero ideal entropy')
# Ideal associated solution A+B <-> AB, total A=B=0.5 mol atoms.
keq=9.;z=(1-sqrt(1/(1+keq)))/2
na=nb=.5-z;nab=z;nt=na+nb+nab
assert isclose((nab/nt)/((na/nt)*(nb/nt)),keq,abs_tol=1e-12)
assert isclose(na+nab,.5) and isclose(nb+nab,.5)
print(f'Association K=9: n_AB={nab:.9f}, species fraction y_AB={nab/nt:.9f}, overall atomic x_B=0.500000000')
# Linearized observations for L(T)=A+B*(T-1000), with known independent Gaussian errors.
observations=[
    ('normalized enthalpy', [1.,-1000.],9100.,100.),
    ('L at 1000K',[1.,0.],5950.,200.),
    ('L at 1500K',[1.,500.],4580.,200.)
]
normal=[[0.,0.],[0.,0.]];rhs=[0.,0.]
for label,row,value,u in observations:
    for i in range(2):
        rhs[i]+=row[i]*value/u**2
        for j in range(2):normal[i][j]+=row[i]*row[j]/u**2
det=normal[0][0]*normal[1][1]-normal[0][1]**2
assert det>0
cov=[[normal[1][1]/det,-normal[0][1]/det],[-normal[1][0]/det,normal[0][0]/det]]
theta=[sum(cov[i][j]*rhs[j] for j in range(2)) for i in range(2)]
chi2=sum((sum(row[i]*theta[i] for i in range(2))-value)**2/u**2 for _,row,value,u in observations)
print(f'Weighted fit: L(1000)={theta[0]:.6f} J/mol, slope={theta[1]:.9f} J/(mol K), chi2={chi2:.9f}')
print(f'Known-noise covariance: u_L1000={sqrt(cov[0][0]):.6f}, u_slope={sqrt(cov[1][1]):.9f}, correlation={cov[0][1]/sqrt(cov[0][0]*cov[1][1]):.9f}')
for temp in (1000.,1500.,3000.):
    row=[1.,temp-1000.]
    pred=sum(row[i]*theta[i] for i in range(2))
    variance=sum(row[i]*cov[i][j]*row[j] for i in range(2) for j in range(2))
    print(f'Conditional model prediction T={temp:.0f} K: L={pred:.6f}, standard uncertainty={sqrt(variance):.6f} J/mol')
# Check the exact weighted gradient, not only a printed optimizer status.
grad=[sum(row[i]*(sum(row[j]*theta[j] for j in range(2))-value)/u**2 for _,row,value,u in observations) for i in range(2)]
assert max(abs(v) for v in grad)<1e-10
# Every binary/ternary boundary misses an interior-only quaternary term.
hidden=-1e6*(.25**4)
assert -1e6*0*.2*.3*.5==0.
print(f'Quaternary-only term: zero on every lower-order boundary, equimolar shift={hidden:.6f} J/mol atoms')
print(f'Five sublattices: 2 constituents each -> {2**5} end members; 4 each -> {4**5}')
print(f'Four components: {comb(4,2)} binaries, {comb(4,3)} ternaries; subsystem counts do not prove extrapolation accuracy')
known={('FE','VA'),('CR','VA'),('FE','C')}
all_members={(a,b) for a in ('FE','CR') for b in ('C','VA')}
missing=sorted(all_members-known)
assert missing==[('CR','C')]
print('Merged (Fe,Cr)1(C,Va)1: missing end member Cr:C must be explicitly assessed or labeled provisional')
print('All model-constraint, weighted-gradient and balance checks PASS; no MCMC or commercial Ag-Cu reproduction claimed')
