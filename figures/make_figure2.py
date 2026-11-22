from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size':7,'axes.linewidth':0.7,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})

ROOT = Path(__file__).resolve().parent
DATA = pd.read_csv(ROOT / "Figure2_source_data.tsv", sep="\t")
FACTOR_ORDER = ["F15", "F9", "F14", "F17", "F13", "F11", "F24", "F20", "F18", "F26", "F23", "F22", "F25"]
GROUP_ORDER = ["lipid_lipoprotein", "fatty_acid", "amino_acid", "energy_ketone", "inflammation_other"]
GROUP_COLORS = {'lipid_lipoprotein':'#7089A6','fatty_acid':'#77A9A2','amino_acid':'#9B8BB4','energy_ketone':'#C69A72','inflammation_other':'#A6A6A6'}
GROUP_LABELS = {'lipid_lipoprotein':'Lipids and lipoproteins','fatty_acid':'Fatty acids','amino_acid':'Amino acids','energy_ketone':'Energy and ketone metabolism','inflammation_other':'Inflammation and other'}
FACTOR_COLORS = {'F9':'#547B9B','F14':'#6B9B8D','F15':'#7388A7','F17':'#9A7FA8'}
CMAP = LinearSegmentedColormap.from_list('muted_div', ['#526F8F','#F6F3EC','#B77474'])

def panel(ax, letter, x=-0.08, y=1.02):
    ax.text(x,y,letter,transform=ax.transAxes,fontsize=9,fontweight='bold',ha='left',va='bottom')

load = DATA[DATA.record_type.eq('loading')].copy()
trait_order = load.drop_duplicates('trait').sort_values(['family','trait'])['trait'].tolist()
# preserve the order encoded in the source file
trait_order = load['trait'].drop_duplicates().tolist()
mat = load.pivot(index='trait',columns='factor',values='display_loading').loc[trait_order,FACTOR_ORDER]
meta = load.drop_duplicates('trait').set_index('trait').loc[trait_order]

fig = plt.figure(figsize=(7.20,7.15), facecolor='white')
gs = fig.add_gridspec(2,2,height_ratios=[2.25,1.0],width_ratios=[1.05,1.55],hspace=0.30,wspace=0.32)
axA = fig.add_subplot(gs[0,:])
im = axA.imshow(mat.values, cmap=CMAP, vmin=-1, vmax=1, aspect='auto', interpolation='nearest')
axA.set_xticks(range(len(FACTOR_ORDER)))
axA.set_xticklabels(FACTOR_ORDER, rotation=0, fontsize=7)
axA.xaxis.tick_top()
axA.tick_params(axis='x',length=0,pad=3)
axA.set_yticks(range(len(trait_order)))
axA.set_yticklabels(trait_order, fontsize=4.8)
axA.tick_params(axis='y',length=0,pad=2)
for sp in axA.spines.values(): sp.set_visible(False)
families = meta['family'].tolist()
start=0
for i in range(1,len(families)+1):
    if i==len(families) or families[i]!=families[start]:
        axA.add_patch(mpl.patches.Rectangle((-1.25,start-0.5),0.18,i-start,facecolor=GROUP_COLORS[families[start]],edgecolor='none',clip_on=False))
        if i<len(families): axA.axhline(i-0.5,color='white',lw=1.2)
        start=i
for x in [2.5,3.5]: axA.axvline(x,color='#3D3D3D',lw=0.55,alpha=0.55)
panel(axA,'A',x=-0.06,y=1.06)
cax=fig.add_axes([0.83,0.922,0.12,0.012])
cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-1,0,1])
cb.set_label('Column-normalized loading',fontsize=6,labelpad=1)
cb.ax.tick_params(labelsize=5,length=2,pad=1)
cb.outline.set_linewidth(0.4)
family_handles=[Patch(facecolor=GROUP_COLORS[g],edgecolor='none',label=GROUP_LABELS[g]) for g in GROUP_ORDER]
axA.legend(handles=family_handles,loc='lower left',bbox_to_anchor=(0.0,1.055),ncol=3,
           fontsize=5.2,handlelength=1.1,handletextpad=.35,columnspacing=.85,borderaxespad=0)

axB=fig.add_subplot(gs[1,0])
met=DATA[DATA.record_type.eq('factor_metric')].set_index('factor').loc[FACTOR_ORDER]
x=met['effective_traits'].astype(float).values
y=met['leading_trait_share'].astype(float).values
single=met['profile_type'].eq('Single-trait').values
axB.scatter(x[~single],y[~single],s=28,color='#728EA3',edgecolor='white',linewidth=0.5,zorder=3)
axB.scatter(x[single],y[single],s=33,color='#B47A7A',marker='D',edgecolor='white',linewidth=0.5,zorder=3)
label_offsets={'F18':(.35,.020),'F20':(.35,-.022),'F24':(.35,.025)}
for fac,xi,yi in zip(FACTOR_ORDER,x,y):
    if fac in {'F22','F23','F25'}:
        target_y={'F22':1.005,'F23':.945,'F25':.885}[fac]
        axB.annotate(fac,xy=(xi,yi),xytext=(3.8,target_y),textcoords='data',fontsize=6.0,
                     va='center',ha='left',color='#333333',
                     arrowprops=dict(arrowstyle='-',color='#A6A6A6',lw=.45,shrinkA=1,shrinkB=2))
        continue
    dx,dy=label_offsets.get(fac,(.35,0.0))
    axB.text(xi+dx,yi+dy,fac,fontsize=6.0,va='center',color='#333333')
axB.set_xlabel('Effective number of contributing traits')
axB.set_ylabel('Leading-trait squared-loading share')
axB.set_xlim(0,27); axB.set_ylim(0,1.04)
axB.set_yticks([0,.25,.5,.75,1.0])
axB.grid(axis='y',color='#E6E6E6',lw=0.5,zorder=0)
axB.scatter([],[],s=28,color='#728EA3',label='Multi-trait')
axB.scatter([],[],s=33,color='#B47A7A',marker='D',label='Single-trait')
axB.legend(loc='upper right',fontsize=5.8,handletextpad=0.4,borderaxespad=0.2)
panel(axB,'B',x=-0.22,y=1.03)

sub=gs[1,1].subgridspec(2,2,hspace=0.55,wspace=0.45)
for i,fac in enumerate(['F9','F14','F15','F17']):
    ax=fig.add_subplot(sub[i//2,i%2])
    d=DATA[(DATA.record_type.eq('representative_loading'))&(DATA.factor.eq(fac))].sort_values('raw_loading')
    yy=np.arange(len(d))
    vals=d['raw_loading'].astype(float).values
    ax.hlines(yy,0,vals,color='#C7CDD2',lw=1)
    ax.scatter(vals,yy,s=22,color=FACTOR_COLORS[fac],edgecolor='white',linewidth=0.4,zorder=3)
    ax.axvline(0,color='#777777',lw=0.5)
    ax.set_yticks(yy); ax.set_yticklabels(d['trait'],fontsize=5.7)
    lim=max(0.3,np.max(np.abs(vals))*1.18); ax.set_xlim(-lim,lim)
    ax.set_xticks([-round(lim*.8,1),0,round(lim*.8,1)])
    ax.tick_params(axis='both',length=2,labelsize=5.2)
    ax.text(0.02,1.05,fac,transform=ax.transAxes,fontsize=6.5,fontweight='bold',color=FACTOR_COLORS[fac])
    if i//2==1: ax.set_xlabel('Loading',fontsize=6)
    else: ax.set_xticklabels([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_color('#8A8A8A')
    if i==0: panel(ax,'C',x=-0.38,y=1.20)

for ext in ['pdf','svg','png']:
    fig.savefig(ROOT/f'Figure2.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight',facecolor='white')
plt.close(fig)
