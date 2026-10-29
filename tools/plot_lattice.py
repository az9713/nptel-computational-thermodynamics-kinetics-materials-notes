"""Original constrained caloric branches and oscillator stiffness effect."""
import os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
os.environ['MPLCONFIGDIR']=str(ROOT/'work/matplotlib-cache');sys.path.insert(0,str(ROOT/'code'))
from ch10_lattice_caloric import solid_hs,liquid_hs,oscillator,R
from math import log
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('dark_background')
fig,axes=plt.subplots(1,2,figsize=(12,4.8),facecolor='#101827')
for a in axes:a.set_facecolor('#101827');a.grid(alpha=.15)
ts=list(range(400,1401,5))
for func,label,color in [(solid_hs,'Solid branch','#7dd3fc'),(liquid_hs,'Liquid branch','#f5a97f')]:
 axes[0].plot(ts,[(func(t)[0]-t*func(t)[1])/1000 for t in ts],label=label,color=color)
axes[0].axvline(1000,color='#a6da95',ls='--',lw=1);axes[0].set(xlabel='Temperature / K',ylabel='Common-reference G / (kJ/mol)',title='Declared branches meet at1000 K');axes[0].legend()
ts=[100*10**(i/80) for i in range(201)]
du=[];df=[]
for t in ts:
 u1,_,f1=oscillator(300,t);u2,_,f2=oscillator(600,t);du.append((u2-u1)/(R*t));df.append((f2-f1)/(R*t))
axes[1].semilogx(ts,du,label='Energy difference / RT',color='#7dd3fc');axes[1].semilogx(ts,df,label='Helmholtz difference / RT',color='#f5a97f');axes[1].axhline(log(2),ls='--',lw=1,color='#a6da95',label='High-temperature ln2 limit for delta-F/RT');axes[1].set(xlabel='Temperature / K',ylabel='Dimensionless difference per mole of modes',title='Doubling frequency: energy and free energy differ');axes[1].legend(fontsize=7)
fig.tight_layout();fig.savefig(ROOT/'assets/lattice-caloric.svg');fig.savefig(ROOT/'work/reader/ch10-lattice.png',dpi=140);plt.close(fig)
print('Wrote original caloric and oscillator curves.')
