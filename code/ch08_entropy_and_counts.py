"""Original Carnot ledger, entropy balances and exact fair-coin counts."""
from math import comb,log,log10,sqrt,isclose
R=8.31446261815324

def coin_probability(n,heads):return comb(n,heads)/2**n

def main():
    hot,cold,ratio=600.,300.,2.
    qhot=R*hot*log(ratio)
    qcold=-R*cold*log(ratio)
    work=qhot+qcold
    assert isclose(qhot/hot+qcold/cold,0,abs_tol=1e-12)
    assert isclose(work/qhot,1-cold/hot)
    # Two finite equal heat-capacity bodies, insulated as a composite.
    capacity=1000.
    final=(hot+cold)/2
    change_hot=capacity*log(final/hot)
    change_cold=capacity*log(final/cold)
    assert change_hot<0<change_cold and change_hot+change_cold>0
    # Fixed-temperature one-step expansion: actual heat differs from reversible heat.
    temperature=300.
    delta_s=R*log(2)
    q_actual=R*temperature/2
    generated=delta_s-q_actual/temperature
    assert generated>0
    n=1000
    middle=coin_probability(n,500)
    band=sum(comb(n,k) for k in range(450,551))/2**n
    assert isclose(sum(coin_probability(4,k) for k in range(5)),1)
    assert middle<.026 and band>.998
    # Two four-oscillator groups share twelve integer energy quanta.
    ways=lambda q:comb(q+3,3)
    weights=[ways(q)*ways(12-q) for q in range(13)]
    assert weights==weights[::-1] and weights.index(max(weights))==6
    assert sum(weights)==comb(12+8-1,8-1)
    print('Carnot 600/300 K, isothermal ratio2: Qhot %.6f J; Qcold %.6f J; work %.6f J; efficiency %.6f.' % (qhot,qcold,work,work/qhot))
    print('Equal finite bodies: final %.1f K; DeltaS_hot %.9f; DeltaS_cold %.9f; total %.9f J/K.' % (final,change_hot,change_cold,change_hot+change_cold))
    print('One-step 300 K expansion: DeltaS_gas %.9f; reservoir %.9f; generated %.9f J/K.' % (delta_s,-q_actual/temperature,generated))
    print('Ideal isolated doubling: Q=0, W=0; DeltaS=R ln2 = %.9f J/K.' % delta_s)
    print('1000 fair coins, log10 multiplicities at 100 / 300 / 500 heads: '+', '.join('%.9f' % log10(comb(n,k)) for k in [100,300,500]))
    print('1000 coins: exact middle probability %.9f; 450..550 band %.9f; fractional standard deviation %.9f.' % (middle,band,1/(2*sqrt(n))))
    print('All-head probability log10 = %.9f; this event has positive probability.' % (-n*log10(2)))
    print('Two four-oscillator groups: modal split6/6; probability %.9f; total microstates %d.' % (weights[6]/sum(weights),sum(weights)))
    print('Checks: cyclic entropy, finite-body signs, actual reservoir heat, exact normalization and multiplicity convolution PASS.')

if __name__=='__main__':main()
