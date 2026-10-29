"""Original Legendre-tangent geometry and declared phase difference."""
import os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
os.environ['MPLCONFIGDIR']=str(ROOT/'work/matplotlib-cache');sys.path.insert(0,str(ROOT/'code'))
from ch09_potentials_entropy import C,TREF,VREF,energy,phase_difference
from math import log
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('dark_background')
fig,axes=plt.subplots(1,2,figsize=(12,4.8),facecolor='#101827')
for ax in axes:ax.set_facecolor('#101827');ax.grid(alpha=.15)
s=[-10+i*.15 for i in range(201)]
axes[0].plot(s,[energy(x,VREF)/1000 for x in s],label='U(S), fixed V',color='#7dd3fc',lw=2.5)
for t,color in [(300,'#a6da95'),(450,'#f5a97f'),(600,'#c6a0f6')]:
 ss=C*log(t/TREF);aa=C*t-t*ss
 axes[0].plot(s,[(t*x+aa)/1000 for x in s],label=f'T={t} K; intercept A={aa/1000:.2f} kJ',color=color,lw=1)
 axes[0].scatter([ss],[C*t/1000],color=color,s=20)
axes[0].set(xlabel='Entropy S / (J/K), chosen reference',ylabel='Energy / kJ',ylim=(-4,15),title='Tangent slope T; vertical intercept U-TS');axes[0].axvline(0,color='white',lw=.6);axes[0].legend(fontsize=7)
temps=list(range(400,1601,4));dg=[phase_difference(t)[2]/1000 for t in temps]
axes[1].plot(temps,dg,color='#7dd3fc',label='Exact declared constant delta-Cp model')
axes[1].plot(temps,[(10000-10*t)/1000 for t in temps],ls='--',color='#f5a97f',label='Freezing delta-H and delta-S at1000 K')
axes[1].axhline(0,color='white',lw=.6);axes[1].scatter([1000],[0],color='#a6da95')
axes[1].set(xlabel='Temperature / K',ylabel='g(beta)-g(alpha) / (kJ/mol)',title='Same crossing; different curvature away from it');axes[1].legend(fontsize=7)
fig.tight_layout();fig.savefig(ROOT/'assets/potentials-entropy.svg');fig.savefig(ROOT/'work/reader/ch09-potentials.png',dpi=140);plt.close(fig)
print('Wrote original potential geometry and heat-capacity difference curves.')
