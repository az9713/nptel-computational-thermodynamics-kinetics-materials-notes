from pathlib import Path
import os
import sys
root=Path(__file__).resolve().parent.parent
(root/'work').mkdir(exist_ok=True)
sys.path.insert(0,str(root))
from numerical_lab import free_energy,coexistence,check_model
os.environ['MPLCONFIGDIR']=str(root/'work/matplotlib-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
x=[i/1000 for i in range(1001)]
y=[free_energy(t,3.0) for t in x]
a,b=coexistence(3.0);g=free_energy(a,3.0)
fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
fig.patch.set_facecolor('#11151d');ax.set_facecolor('#11151d')
ax.plot(x,y,color='#78c8f7',label='Homogeneous regular solution')
ax.hlines(g,a,b,color='#f2bf78',linestyle='--',label='Two-phase common tangent')
ax.scatter([a,b],[g,g],color='#f2bf78',zorder=3)
ax.set(xlabel='Component B mole fraction x',ylabel='Mixing Gibbs energy / (RT)',xlim=(0,1))
ax.tick_params(colors='#dce4ef');ax.xaxis.label.set_color('#dce4ef');ax.yaxis.label.set_color('#dce4ef')
for s in ax.spines.values():s.set_color('#7f8997')
ax.grid(alpha=.15);ax.legend(facecolor='#202733',labelcolor='#dce4ef')
fig.savefig(root/'assets/regular-solution.svg',metadata={'Title':'Original synthetic regular-solution construction'})
fig.savefig(root/'work/regular-solution.png',dpi=120)
print(check_model())
