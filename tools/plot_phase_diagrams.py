"""Plot synthetic equilibrium diagram from the chapter's declared model."""
import os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
os.environ['MPLCONFIGDIR']=str(ROOT/'work/matplotlib-cache')
sys.path.insert(0,str(ROOT/'code'))
from ch03_phase_diagrams import endpoints,invert_endpoint,TMA,TMB
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('dark_background')
T=[TMA+i*(TMB-TMA)/400 for i in range(401)]
xl,xs=zip(*(endpoints(t) for t in T))
x0=.45;ts=invert_endpoint(x0,1);tl=invert_endpoint(x0,0)
fig,(ax,bx)=plt.subplots(1,2,figsize=(11,5.5),gridspec_kw={'width_ratios':[1.5,1]},facecolor='#101827')
for a in (ax,bx):a.set_facecolor('#101827');a.grid(alpha=.16)
ax.plot(xl,T,color='#f5a97f',label='Liquidus (liquid endpoint)')
ax.plot(xs,T,color='#7dd3fc',label='Solidus (solid endpoint)')
ax.fill_betweenx(T,xl,xs,color='#c4b5fd',alpha=.14)
ax.axvline(x0,color='#c4b5fd',ls='--',lw=1)
a,b=endpoints(1200)
ax.plot([a,b],[1200,1200],color='white',lw=2)
ax.scatter([a,b],[1200,1200],color=['#f5a97f','#7dd3fc'])
ax.text(.12,1360,'Liquid',fontsize=12);ax.text(.72,1150,'Solid',fontsize=12)
ax.annotate('Tie line at 1200 K',xy=((a+b)/2,1200),xytext=(.06,1070),arrowprops={'arrowstyle':'->','color':'white'})
ax.set(xlabel='Atomic fraction B',ylabel='Temperature (K)',xlim=(0,1),ylim=(950,1550),title='Synthetic ideal solid/liquid model')
ax.legend(loc='upper left',fontsize=9)
fractions=[]
for t in T:
 if t<=ts:f=1
 elif t>=tl:f=0
 else:
  a,b=endpoints(t);f=(x0-a)/(b-a)
 fractions.append(f)
bx.plot(fractions,T,color='#c4b5fd',lw=2)
bx.axhline(ts,color='#7dd3fc',ls=':');bx.axhline(tl,color='#f5a97f',ls=':')
bx.set(xlabel='Solid atomic-amount fraction',ylabel='Temperature (K)',xlim=(-.05,1.05),ylim=(1100,1300),title='Fixed overall x_B = 0.45')
fig.suptitle('Equilibrium amounts are functions of state, not a cooling-time history',fontsize=13)
fig.tight_layout();fig.savefig(ROOT/'assets/isomorphous-model.svg');fig.savefig(ROOT/'work/reader/ch03-model.png',dpi=130);plt.close(fig)
print('Wrote declared-model diagram and equilibrium amount curve.')
