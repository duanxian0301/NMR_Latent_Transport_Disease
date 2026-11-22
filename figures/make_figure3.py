from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans','Liberation Sans']; plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'font.size':7,'axes.linewidth':0.7,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
ROOT=Path(__file__).resolve().parent
D=pd.read_csv(ROOT/'Figure3_source_data.tsv',sep='\t')
GROUP_COLORS={'lipid_lipoprotein':'#7089A6','fatty_acid':'#77A9A2','amino_acid':'#9B8BB4','energy_ketone':'#C69A72','inflammation_other':'#A6A6A6'}
GROUP_LABELS={'lipid_lipoprotein':'Lipids and lipoproteins','fatty_acid':'Fatty acids','amino_acid':'Amino acids','energy_ketone':'Energy and ketone metabolism','inflammation_other':'Inflammation and other'}
CMAP=LinearSegmentedColormap.from_list('sim',['#F4F2EC','#AFC2CF','#526F8F'])
FACTORS=[f'F{i}' for i in range(1,31)]
KARJ=[f'Karj_F{i}' for i in range(1,16)]

def panel(ax,l,x=-0.12,y=1.02): ax.text(x,y,l,transform=ax.transAxes,fontsize=9,fontweight='bold',ha='left',va='bottom')

fig=plt.figure(figsize=(7.2,7.2),facecolor='white')
gs=fig.add_gridspec(2,2,width_ratios=[1.0,1.45],height_ratios=[1.2,1.0],hspace=0.34,wspace=0.34)

axA=fig.add_subplot(gs[0,0])
a=D[D.record_type.eq('trait_rg')].sort_values('estimate')
y=np.arange(len(a)); est=a.estimate.astype(float).values; lo=a.ci_low.astype(float).values; hi=a.ci_high.astype(float).values
for yi,(l,h,c,u) in enumerate(zip(lo,hi,a.family,a.usable)):
    col=GROUP_COLORS.get(c,'#888888') if str(u).lower()=='true' else '#B55D5D'
    axA.plot([l,h],[yi,yi],color=col,lw=0.75,alpha=.8)
    axA.scatter(est[yi],yi,s=8,color=col,zorder=3,marker='o' if str(u).lower()=='true' else 'x',linewidth=.7)
axA.axvline(0,color='#777777',lw=.6,ls='--')
axA.set_yticks([]); axA.set_xlabel('Cross-resource genetic correlation')
axA.set_xlim(min(-.05,np.nanmin(lo)-.04),max(1.10,np.nanmax(hi)+.02))
axA.grid(axis='x',color='#E9E9E9',lw=.5)
axA.text(.02,.03,'61 metabolites; 60 usable',transform=axA.transAxes,fontsize=5.8,color='#555555')
for g,col in GROUP_COLORS.items(): axA.scatter([],[],s=12,color=col,label=GROUP_LABELS[g])
axA.legend(loc='lower left',bbox_to_anchor=(0,1.01),ncol=2,fontsize=5.2,handletextpad=.3,columnspacing=.7)
panel(axA,'A',x=-.18,y=1.13)

axB=fig.add_subplot(gs[0,1])
b=D[D.record_type.eq('factor_trait_rg')].copy()
order=['F9','F11','F13','F14','F15','F17','F18','F20','F22','F23','F24','F25','F26']
for j,fac in enumerate(order[::-1]):
    dd=b[b.factor.eq(fac)]
    for entity,off,col,mark in [('Projected factor',-.12,'#567E98','o'),('Dominant metabolite',.12,'#9D7A9B','D')]:
        r=dd[dd.entity.eq(entity)].iloc[0]; yy=j+off
        est=float(r.estimate); lo=float(r.ci_low); hi=float(r.ci_high); usable=str(r.usable).lower()=='true'
        axB.plot([lo,hi],[yy,yy],color=col if usable else '#B9B9B9',lw=1.0,zorder=1)
        axB.scatter(est,yy,s=21,color=col if usable else 'white',edgecolor=col if usable else '#777777',marker=mark,linewidth=.7,zorder=2)
        if not usable: axB.scatter(est,yy,s=11,color='#777777',marker='x',linewidth=.7,zorder=3)
axB.axvline(0,color='#777777',lw=.6,ls='--')
axB.set_yticks(range(len(order))); axB.set_yticklabels(order[::-1],fontsize=6)
axB.set_xlabel('Cross-resource genetic correlation')
axB.set_xlim(-.35,1.55); axB.grid(axis='x',color='#E9E9E9',lw=.5)
axB.scatter([],[],s=21,color='#567E98',label='Projected factor')
axB.scatter([],[],s=21,color='#9D7A9B',marker='D',label='Dominant metabolite')
axB.scatter([],[],s=18,color='white',edgecolor='#777777',label='Not usable')
axB.legend(loc='lower left',bbox_to_anchor=(1.01,0.0),fontsize=5.8,ncol=1,borderaxespad=0)
panel(axB,'B',x=-.12,y=1.03)

axC=fig.add_subplot(gs[1,0])
c=D[D.record_type.eq('loading_similarity')]
M=c.pivot(index='karj_factor',columns='factor',values='estimate').reindex(index=KARJ,columns=FACTORS)
im=axC.imshow(M.values,cmap=CMAP,vmin=0,vmax=1,aspect='auto',interpolation='nearest')
axC.set_xticks(range(30)); axC.set_xticklabels(FACTORS,rotation=90,fontsize=4.9)
axC.set_yticks(range(15)); axC.set_yticklabels([x.replace('Karj_','K') for x in KARJ],fontsize=5.3)
axC.tick_params(length=0)
for _,r in c[c.is_unique_match.astype(str).str.lower().eq('true')].iterrows():
    x=FACTORS.index(r.factor); y0=KARJ.index(r.karj_factor)
    ec='#9A5F77' if float(r.estimate)>=.85 else '#5F5F5F'
    axC.add_patch(mpl.patches.Rectangle((x-.47,y0-.47),.94,.94,fill=False,edgecolor=ec,lw=.75))
for sp in axC.spines.values(): sp.set_visible(False)
panel(axC,'C',x=-.17,y=1.03)
cax=axC.inset_axes([.56,1.035,.38,.035]); cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[0,.5,1]); cb.set_label('Absolute loading cosine',fontsize=5.3,labelpad=1); cb.ax.tick_params(labelsize=4.8,length=1.5,pad=1); cb.outline.set_linewidth(.35)

subD=gs[1,1].subgridspec(2,1,height_ratios=[1.0,.78],hspace=.48)
axD=fig.add_subplot(subD[0,0])
d=D[D.record_type.eq('subspace_reference')]
ypos=[1.35,.75]
for yi,(_,r) in zip(ypos,d.iterrows()):
    axD.plot([float(r.null_mean),float(r.estimate)],[yi,yi],color='#C3C3C3',lw=2,zorder=1)
    axD.scatter(float(r.null_mean),yi,s=35,color='#A9A9A9',marker='s',zorder=2)
    axD.scatter(float(r.estimate),yi,s=42,color='#587B98',zorder=3)
    axD.text(max(float(r.null_mean),float(r.estimate))+.012,yi,f"P = {float(r.p_value):.3f}",va='center',fontsize=5.8)
axD.set_yticks(ypos); axD.set_yticklabels(['UKB residual reference','Karjalainen residual reference'],fontsize=6)
axD.set_xlabel('Mean canonical correlation')
axD.set_xlim(.81,.925); axD.set_ylim(.48,1.65); axD.grid(axis='x',color='#E9E9E9',lw=.5)
axD.scatter([],[],s=38,color='#587B98',label='Observed')
axD.scatter([],[],s=32,color='#A9A9A9',marker='s',label='Null mean')
axD.legend(loc='upper left',fontsize=5.8,ncol=2)
panel(axD,'D',x=-.13,y=1.03)

ax2=fig.add_subplot(subD[1,0])
low=D[D.record_type.eq('low_dimension')].sort_values('k')
ax2.plot(low.k.astype(int),low.estimate.astype(float),color='#8F789B',lw=1.2,marker='o',ms=3.5)
ax2.set_xticks(low.k.astype(int)); ax2.set_xlabel('Leading dimensions',fontsize=6,labelpad=1); ax2.set_ylabel('Mean canonical correlation',fontsize=6,labelpad=2)
ax2.tick_params(labelsize=5.3,length=2); ax2.set_ylim(.38,.66); ax2.grid(axis='y',color='#E9E9E9',lw=.4)

for ext in ['pdf','svg','png']:
    fig.savefig(ROOT/f'Figure3.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight',facecolor='white')
plt.close(fig)
