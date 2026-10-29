from pathlib import Path
import os,runpy,io,contextlib
ROOT=Path(__file__).resolve().parent.parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib-cache'))
with contextlib.redirect_stdout(io.StringIO()):v=runpy.run_path(str(ROOT/'code/ch19_diffusion_mobility.py'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'text.color':'#e7edf5','axes.labelcolor':'#e7edf5','xtick.color':'#c3d0df','ytick.color':'#c3d0df','axes.edgecolor':'#8394a8','font.size':11,'svg.fonttype':'none','svg.hashsalt':'thermo-ch19'})
fig,axs=plt.subplots(1,2,figsize=(12,5.2),facecolor='#111827')
for ax in axs:ax.set_facecolor('#111827');ax.grid(alpha=.15,color='#94a3b8')
x=np.linspace(0,.003,401)
for factor,color in [(.25,'#fbbf24'),(1,'#67e8f9'),(4,'#c4b5fd')]:
 t=v['time']*factor
 axs[0].plot(x*1000,[v['concentration'](xx,t) for xx in x],color=color,label=f'{t/3600:.2f} h')
axs[0].scatter([1],[3000],color='#e7edf5',s=25)
axs[0].set(xlabel='Depth (mm)',ylabel='Concentration (mol/m³)',title='Maintained surface reservoir; constant D')
axs[0].legend(facecolor='#182335',labelcolor='#e7edf5')
ratio=np.linspace(0,1.65,400);k=ratio*v['kstar'];sigma=-v['mob']*k*k*(v['h']+v['kap']*k*k)
axs[1].plot(ratio,sigma*1e4,color='#67e8f9')
axs[1].axhline(0,color='#c3d0df',lw=.7);axs[1].axvline(1,color='#fbbf24',ls='--',lw=.9)
axs[1].set(xlabel='Wavenumber / fastest-growing wavenumber',ylabel='Linear growth rate (10⁻⁴/s)',title=f"Fastest wavelength {v['wavelength']*1e9:.2f} nm")
fig.tight_layout();fig.savefig(ROOT/'assets/diffusion-mobility.svg',metadata={'Date':None})
if __name__=='__main__':print('Wrote assets/diffusion-mobility.svg')
