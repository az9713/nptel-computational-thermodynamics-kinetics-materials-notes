"""Original regular-solution mixing curve and consistent activities."""
import os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
os.environ['MPLCONFIGDIR']=str(ROOT/'work/matplotlib-cache');sys.path.insert(0,str(ROOT/'code'))
from ch11_mixtures_activity import molar_g,chemical,GA,GB,R
from math import exp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('dark_background');fig,axes=plt.subplots(1,2,figsize=(12,4.8),facecolor='#101827')
for a in axes:a.set_facecolor('#101827');a.grid(alpha=.15)
xs=[.001+i*.998/500 for i in range(501)]
for t,c in [(600,'#f5a97f'),(800,'#7dd3fc')]:
 axes[0].plot(xs,[(molar_g(x,t)-(1-x)*GA-x*GB)/1000 for x in xs],label=f'T={t} K',color=c)
axes[0].axhline(0,color='white',lw=.6);axes[0].set(xlabel='Mole fraction B',ylabel='Mixing Gibbs energy / (kJ/mol)',title='Negative mixing energy is not global stability');axes[0].legend()
for index,label,c in [(0,'A','#a6da95'),(1,'B','#c6a0f6')]:
 act=[exp((chemical(x,800)[index]-[GA,GB][index])/(R*800)) for x in xs]
 axes[1].plot(xs,act,color=c,label=f'Regular a{label}')
 axes[1].plot(xs,[1-x if index==0 else x for x in xs],ls='--',color=c,label=f'Ideal a{label}')
axes[1].set(xlabel='Mole fraction B',ylabel='Dimensionless matched-reference activity',title='Partial derivatives generate both activities');axes[1].legend(fontsize=8)
fig.tight_layout();fig.savefig(ROOT/'assets/mixtures-activity.svg');fig.savefig(ROOT/'work/reader/ch11-mixtures.png',dpi=140);plt.close(fig)
print('Wrote original mixing-energy and activity curves.')
