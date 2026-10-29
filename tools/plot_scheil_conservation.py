from pathlib import Path
import os
ROOT=Path(__file__).resolve().parent.parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
c0=.01;k=.13;ce=.12;rE=(ce/c0)**(1/(k-1));fE=1-rE
f=np.linspace(.0001,fE,600);r=1-f;cl=c0*r**(k-1);cs=k*cl;avg=c0*(1-r**k)/f
plt.rcParams.update({'text.color':'#e7edf5','axes.labelcolor':'#e7edf5','xtick.color':'#c3d0df','ytick.color':'#c3d0df','axes.edgecolor':'#8394a8','font.size':11,'svg.fonttype':'none','svg.hashsalt':'thermo-ch16'})
fig,axs=plt.subplots(1,2,figsize=(12,5.2),facecolor='#111827')
for ax in axs:ax.set_facecolor('#111827');ax.grid(alpha=.15,color='#94a3b8')
ax=axs[0];ax.plot(f,100*cl,label='Remaining liquid',color='#67e8f9');ax.plot(f,100*cs,label='New solid',color='#fbbf24');ax.plot(f,100*avg,label='Mean frozen solid',color='#f472b6');ax.axvline(fE,color='#c3d0df',ls='--',alpha=.8);ax.set(xlabel='Solid mass fraction before eutectic',ylabel='Solute mass percent',title='Different meanings of solid composition',xlim=(0,1));ax.legend(facecolor='#182335',labelcolor='#e7edf5',fontsize=9)
ax=axs[1];ff=np.linspace(0,1,600);eqcl=c0/(k+(1-k)*(1-ff));ax.plot(ff,933-600*eqcl,label='Equilibrium: homogeneous solid',color='#fbbf24');ax.plot(f,933-600*cl,label='Scheil primary-solid path',color='#67e8f9');ax.plot([fE,1],[861,861],label='Terminal two-solid mixture',color='#f472b6',lw=3);ax.scatter([fE],[861],color='#e7edf5');ax.set(xlabel='Total solid mass fraction',ylabel='Temperature (K)',title='A path without a cooling clock',xlim=(0,1),ylim=(850,935));ax.legend(facecolor='#182335',labelcolor='#e7edf5',fontsize=9)
fig.suptitle('Original synthetic constant-k model; not an assessed Al–Si calculation',fontsize=13)
fig.tight_layout();fig.savefig(ROOT/'assets/scheil-conservation.svg',metadata={'Date':None})
if __name__=='__main__':print('Wrote assets/scheil-conservation.svg')
