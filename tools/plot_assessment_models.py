from pathlib import Path
import os,runpy,io,contextlib
ROOT=Path(__file__).resolve().parent.parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib-cache'))
with contextlib.redirect_stdout(io.StringIO()):values=runpy.run_path(str(ROOT/'code/ch18_models_and_assessment.py'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'text.color':'#e7edf5','axes.labelcolor':'#e7edf5','xtick.color':'#c3d0df','ytick.color':'#c3d0df','axes.edgecolor':'#8394a8','font.size':11,'svg.fonttype':'none','svg.hashsalt':'thermo-ch18'})
fig,axs=plt.subplots(1,2,figsize=(12,5.2),facecolor='#111827')
for ax in axs:ax.set_facecolor('#111827');ax.grid(alpha=.15,color='#94a3b8')
eta=np.linspace(-.999,.999,700);u=(1+eta)/2;q=u*np.log(u)+(1-u)*np.log(1-u)
for ratio,color in [(.8,'#67e8f9'),(1.2,'#fbbf24')]:
 dg=-values['omega']*eta**2/2+values['R']*(ratio*values['tc'])*(q+np.log(2))
 axs[0].plot(eta,dg,label=f'T/Tc={ratio}',color=color)
axs[0].axhline(0,color='#c3d0df',lw=.6);axs[0].set(xlabel='Order parameter eta at x_B=0.5',ylabel='G(eta)-G(0) (J/mol atoms)',title='Internal state minimization');axs[0].legend(facecolor='#182335',labelcolor='#e7edf5')
ts=np.linspace(600,3200,500);vv=np.column_stack((np.ones_like(ts),ts-1000));theta=np.array(values['theta']);cov=np.array(values['cov']);pred=vv@theta;u=np.sqrt(np.einsum('ij,jk,ik->i',vv,cov,vv))
ax=axs[1];ax.plot(ts,pred,color='#67e8f9',label='Weighted linear interaction');ax.fill_between(ts,pred-u,pred+u,color='#67e8f9',alpha=.2,label='Conditional +/-1 standard deviation');ax.scatter([1000,1500],[5950,4580],color='#fbbf24',label='Direct L observations');ax.axvline(1500,color='#c3d0df',ls='--',lw=.8);ax.set(xlabel='Temperature (K)',ylabel='L(T) (J/mol atoms)',title='Parameter-induced uncertainty in a chosen model');ax.legend(facecolor='#182335',labelcolor='#e7edf5',fontsize=9)
fig.tight_layout();fig.savefig(ROOT/'assets/assessment-models.svg',metadata={'Date':None})
if __name__=='__main__':print('Wrote assets/assessment-models.svg')
