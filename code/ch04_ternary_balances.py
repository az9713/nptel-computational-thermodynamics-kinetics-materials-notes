"""Exact barycentric balances on declared atom or mass bases."""
from fractions import Fraction as F
from math import sqrt


def solve3(columns, rhs):
    rows=[[F(columns[j][i]) for j in range(3)]+[F(rhs[i])] for i in range(3)]
    for j in range(3):
        pivot=next((i for i in range(j,3) if rows[i][j]),None)
        if pivot is None:
            raise ValueError('Degenerate tie triangle: amounts not uniquely determined')
        rows[j],rows[pivot]=rows[pivot],rows[j]
        scale=rows[j][j]
        rows[j]=[v/scale for v in rows[j]]
        for i in range(3):
            if i!=j:
                scale=rows[i][j]
                rows[i]=[a-scale*b for a,b in zip(rows[i],rows[j])]
    return tuple(rows[i][3] for i in range(3))


def weighted(columns, fractions):
    return tuple(sum(f*v[i] for v,f in zip(columns,fractions)) for i in range(3))


def cartesian(x):
    return float(x[1]+x[2]/2),sqrt(3)*float(x[2])/2


def twice_area(p,q,r):
    return (q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0])


if __name__=='__main__':
    phases=[(F(4,5),F(3,20),F(1,20)),
            (F(1,10),F(4,5),F(1,10)),
            (F(3,20),F(3,20),F(7,10))]
    amounts=(F(1,5),F(3,10),F(1,2))
    overall=weighted(phases,amounts)
    assert overall==(F(53,200),F(69,200),F(39,100))
    assert solve3(phases,overall)==amounts
    p,q,r=map(cartesian,phases);x=cartesian(overall)
    total=twice_area(p,q,r)
    area_amounts=[twice_area(x,q,r)/total,twice_area(p,x,r)/total,
                  twice_area(p,q,x)/total]
    assert all(abs(a-float(f))<1e-14 for a,f in zip(area_amounts,amounts))
    checks=0
    for i in range(21):
        for j in range(21-i):
            f=(F(i,20),F(j,20),F(20-i-j,20))
            assert solve3(phases,weighted(phases,f))==f
            checks+=1
    outside=solve3(phases,(1,0,0))
    assert min(outside)<0
    print('Synthetic overall (A,B,C):',*[str(v) for v in overall])
    print('Recovered phase atomic-amount fractions:',*[str(v) for v in amounts])
    print('Area ratios:',*[f'{v:.6f}' for v in area_amounts])
    print('Exact grid reconstructions:',checks)
    # A vertical section x_B=.4 contains the bulk but neither tie endpoint.
    ends=[(F(7,10),F(1,10),F(1,5)),(F(1,10),F(3,5),F(3,10))]
    bulk=weighted(ends,(F(2,5),F(3,5)))
    assert bulk==(F(17,50),F(2,5),F(13,50))
    assert all(v[1]!=bulk[1] for v in ends)
    print('Off-section tie-line bulk:',*[str(v) for v in bulk])
    # NIST legacy calculated eutectic table, Ag/Cu/Sn MASS fractions.
    solids=[(F(7317,10000),F(0),F(2683,10000)),
            (F(0),F(3907,10000),F(6093,10000)),
            (F(7,10000),F(0),F(9993,10000))]
    liquid=(F(373,10000),F(85,10000),F(9542,10000))
    products=solve3(solids,liquid)
    assert all(f>0 for f in products) and sum(products)==1
    assert weighted(solids,products)==liquid
    print('NIST table all-solid product MASS fractions (Ag3Sn,Cu6Sn5,Sn):',
          *[f'{float(v):.9f}' for v in products])
    # Ternary term invisible on every binary boundary.
    center=(F(1,3),)*3
    for edge in [(F(1,4),F(3,4),F(0)),(F(0),F(1,4),F(3,4)),
                 (F(3,4),F(0),F(1,4))]:
        assert F(2700)*edge[0]*edge[1]*edge[2]==0
    assert F(2700)*center[0]*center[1]*center[2]==100
    print('Center ternary energy for lambda=+/-2700 J/mol: +/-100 J/mol')
    print('All conservation, area, section and boundary checks passed.')