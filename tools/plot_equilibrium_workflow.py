from pathlib import Path
import os
ROOT=Path(__file__).resolve().parent.parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
R=8.3145;T=600.;z=.4
x=np.linspace(.0001,.9999,900);q=(1-x)*np.log(1-x)+x*np.log(x)
a=1/(1+np.exp(6000/(R*T)));b=1-a;tie=R*T*np.log(1-a)
plt.rcParams.update({'text.color':'#e7edf5','axes.labelcolor':'#e7edf5','xtick.color':'#c3d0df','ytick.color':'#c3d0df','axes.edgecolor':'#8394a8','font.size':11,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,2,figsize=(12,5.0),facecolor='#111827')
for ax in axes:ax.set_facecolor('#111827');ax.grid(alpha=.15,color='#94a3b8')
ax=axes[0];ax.plot(x,6000*x+R*T*q,label='Homogeneous ALPHA',color='#67e8f9');ax.plot(x,6000*(1-x)+R*T*q,label='Homogeneous BETA',color='#fbbf24');ax.plot([a,b],[tie,tie],label='Equilibrium tie line',color='#f472b6',lw=3);ax.scatter([a,b],[tie,tie],color='#f472b6');ax.set(xlabel='x(NI), teaching variable',ylabel='Molar Gibbs energy (J/mol)',title='Restricted phase set at 600 K');ax.legend(facecolor='#182335',labelcolor='#e7edf5',fontsize=9)
temps=np.linspace(400,2000,600);aa=1/(1+np.exp(6000/(R*temps)));fb=np.clip((z-aa)/(1-2*aa),0,1)
ax=axes[1];ax.plot(temps,fb,color='#fbbf24',label='Analytic beta fraction');ax.scatter([400,600,1000],[.360586,.314132,.210926],color='#67e8f9',s=35,label='Executed engine values');ax.axvline(6000/(R*np.log(1.5)),color='#a3b4c8',ls='--');ax.set(xlabel='Temperature (K)',ylabel='Beta atomic-mole fraction',title='Fixed bulk x(NI) = 0.4');ax.legend(facecolor='#182335',labelcolor='#e7edf5',fontsize=9);ax.set_xlim(350,2050)
fig.suptitle('Original synthetic model: GAMMA excluded in both panels',color='#e7edf5');fig.tight_layout(rect=(0,.05,1,.95));fig.text(.5,.015,'Analytic curves, with three actually executed equilibrium checks; not a measured Cu-Ni system',ha='center',color='#c3d0df',fontsize=10)
fig.savefig(ROOT/'assets/equilibrium-workflow.svg',facecolor=fig.get_facecolor());fig.savefig(ROOT/'work/author/ch14-figure.png',facecolor=fig.get_facecolor(),dpi=130)
