"""Original supporting-line and regular-solution stability constructions."""
import os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
os.environ['MPLCONFIGDIR']=str(ROOT/'work/matplotlib-cache');sys.path.insert(0,str(ROOT/'code'))
from ch12_global_tangents import ga,gb,g_regular,bisect,R
from math import sqrt,log
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('dark_background');fig,axes=plt.subplots(1,2,figsize=(12,4.8),facecolor='#101827')
for ax in axes:ax.set_facecolor('#101827');ax.grid(alpha=.15)
xs=[i/500 for i in range(501)]
axes[0].plot(xs,[ga(x) for x in xs],label='alpha',color='#7dd3fc');axes[0].plot(xs,[gb(x) for x in xs],label='beta',color='#f5a97f');axes[0].axhline(0,color='#a6da95',lw=2,label='Global supporting line')
axes[0].plot(xs,[2000*x-500 for x in xs],ls='--',color='#c6a0f6',label='Parent tangent at x=0.3');axes[0].plot(xs,[2000*x-1700 for x in xs],ls=':',color='#c6a0f6',label='Parallel product tangent at x=0.9')
axes[0].scatter([.2,.8],[0,0],color='#a6da95');axes[0].scatter([.3],[ga(.3)],color='#7dd3fc');axes[0].set(xlabel='Mole fraction B',ylabel='Declared relative g / (J/mol)',ylim=(-1000,2500),title='Finite equilibrium drop and infant-phase driving differ');axes[0].legend(fontsize=7)
t=600;omega=12000
b=bisect(lambda x:R*t*log(x/(1-x))+omega*(1-2*x),1e-12,.499999)
s=(1-sqrt(1-2*R*t/omega))/2
xi=[.001+i*.998/500 for i in range(501)]
axes[1].plot(xi,[g_regular(x,t) for x in xi],color='#7dd3fc',label='Single regular-solution branch');axes[1].axhline(g_regular(b,t),color='#a6da95',label='Global two-phase line')
for v in [b,1-b]:axes[1].axvline(v,color='#a6da95',ls='--',lw=1)
for v in [s,1-s]:axes[1].axvline(v,color='#f5a97f',ls=':',lw=1)
axes[1].axvspan(s,1-s,color='#f5a97f',alpha=.12,label='Negative-curvature interval')
axes[1].set(xlabel='Mole fraction B',ylabel='Mixing energy / (J/mol)',title='600 K: binodal endpoints enclose spinodal endpoints');axes[1].legend(fontsize=7)
fig.tight_layout();fig.savefig(ROOT/'assets/global-tangents.svg');fig.savefig(ROOT/'work/reader/ch12-tangents.png',dpi=140);plt.close(fig)
print('Wrote original global support and local/global stability comparisons.')
