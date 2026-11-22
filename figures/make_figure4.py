from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans','Liberation Sans']; plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'font.size':7,'axes.linewidth':0.7,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
ROOT=Path(__file__).resolve().parent
D=pd.read_csv(ROOT/'Figure4_source_data.tsv',sep='\t')
FACTOR_ORDER=['F9','F11','F13','F14','F15','F17','F18','F20','F22','F23','F24','F25','F26']
COL={'Detected by at least one metabolite':'#7089A6','Dominant-metabolite threshold shadow':'#C7A67A','Multivariate-threshold candidate':'#9B7D9E'}
FACCOL={'F9':'#587E9B','F15':'#6B9B8D','F17':'#9A7FA8'}
CMAP=LinearSegmentedColormap.from_list('rho',['#637F9D','#F5F2EB','#B87373'])

def panel(ax,l,x=-.10,y=1.02): ax.text(x,y,l,transform=ax.transAxes,fontsize=9,fontweight='bold',ha='left',va='bottom')

fig=plt.figure(figsize=(7.2,6.85),facecolor='white')
gs=fig.add_gridspec(2,2,width_ratios=[1.05,1.45],height_ratios=[.9,1.25],hspace=.36,wspace=.34)

axA=fig.add_subplot(gs[0,0]); a=D[D.record_type.eq('locus_class')]
tot=a.groupby('category')['count'].sum(); total=int(tot.sum()); captured=int(tot['Detected by at least one metabolite']); rem=total-captured
axA.barh([.82],[captured],height=.25,color=COL['Detected by at least one metabolite'],edgecolor='none')
axA.barh([.82],[rem],left=[captured],height=.25,color='#D8D2C8',edgecolor='none')
axA.text(captured/2,.82,f'{captured:,}',ha='center',va='center',color='white',fontsize=7,fontweight='bold')
axA.text(total+12,.82,f'{total:,} factor–lead associations',ha='left',va='center',fontsize=6.2)
shadow=int(tot['Dominant-metabolite threshold shadow']); cand=int(tot['Multivariate-threshold candidate'])
ins=axA.inset_axes([0.0,.08,.72,.28])
ins.barh([0],[shadow],height=.5,color=COL['Dominant-metabolite threshold shadow'],edgecolor='none')
ins.barh([0],[cand],left=[shadow],height=.5,color=COL['Multivariate-threshold candidate'],edgecolor='none')
ins.text(shadow/2,0,str(shadow),ha='center',va='center',fontsize=6.5)
ins.text(shadow+cand/2,0,str(cand),ha='center',va='center',fontsize=6.5,color='white')
ins.set_xlim(0,13); ins.set_ylim(-.5,.5); ins.set_xticks([]); ins.set_yticks([])
for sp in ins.spines.values(): sp.set_visible(False)
axA.text(0,.02,'Magnified: 13 associations not captured by a constituent metabolite',fontsize=5.3,color='#555555',transform=axA.transAxes)
axA.set_xlim(0,total*1.23); axA.set_ylim(0,1.18); axA.set_yticks([]); axA.set_xticks([])
for sp in axA.spines.values(): sp.set_visible(False)
handles=[Patch(facecolor=COL[k],edgecolor='none',label=k) for k in COL]
axA.legend(handles=handles,loc='upper left',bbox_to_anchor=(0,1.05),fontsize=5.4,handlelength=1.1)
panel(axA,'A',x=-.08,y=1.04)

axB=fig.add_subplot(gs[0,1]); b=D[D.record_type.eq('locus_contribution')].copy()
traits=[]
for t in b.sort_values(['factor','rank']).trait:
    if t not in traits: traits.append(t)
rows=[('F9','rs10808546'),('F15','rs2980869'),('F17','rs11718898')]
for yi,(fac,snp) in enumerate(rows[::-1]):
    dd=b[(b.factor.eq(fac))&(b.lead_snp.eq(snp))]
    for _,r in dd.iterrows():
        x=traits.index(r.trait); size=20+145*float(r.relative_to_row_max)
        axB.scatter(x,yi,s=size,color=FACCOL[fac],alpha=.85,edgecolor='white',linewidth=.5)
axB.set_xticks(range(len(traits))); axB.set_xticklabels(traits,rotation=55,ha='right',fontsize=5.2)
axB.set_yticks(range(3)); axB.set_yticklabels([f'{f}  {s}' for f,s in rows[::-1]],fontsize=6)
axB.set_xlim(-.7,len(traits)-.3); axB.set_ylim(-.6,2.6); axB.grid(axis='x',color='#ECECEC',lw=.4)
axB.tick_params(axis='both',length=0)
for sp in axB.spines.values(): sp.set_visible(False)
for s,label in [(45,'0.2'),(95,'0.5'),(165,'1.0')]: axB.scatter([],[],s=s,color='#8A9BA8',alpha=.8,edgecolor='white',label=label)
axB.legend(title='Relative contribution',loc='upper right',fontsize=5,title_fontsize=5.4,ncol=3,handletextpad=.2,columnspacing=.6)
panel(axB,'B',x=-.10,y=1.04)

axC=fig.add_subplot(gs[1,0]); c=D[D.record_type.eq('same_factor_lava')].set_index('factor').loc[FACTOR_ORDER]
y=np.arange(len(c))[::-1]; x=c.n_fdr_within.astype(float).values; sizes=18+1.1*c.n_tests_non_mhc.astype(float).values; vals=c.median_rho_fdr.astype(float).values
sc=axC.scatter(x,y,s=sizes,c=vals,cmap=CMAP,norm=Normalize(-1,1),edgecolor='white',linewidth=.6)
for xi,yi,neg in zip(x,y,c.n_negative_fdr.astype(int)):
    if neg>0: axC.text(xi+7,yi,f'{neg}−',fontsize=5.2,va='center',color='#9A5B5B')
axC.set_yticks(y); axC.set_yticklabels(c.index,fontsize=6)
axC.set_xlabel('Within-factor FDR-significant local correlations')
axC.set_xlim(-5,max(x)+45); axC.grid(axis='x',color='#E9E9E9',lw=.5)
cax=axC.inset_axes([.48,1.035,.48,.032]); cb=fig.colorbar(sc,cax=cax,orientation='horizontal',ticks=[-1,0,1]); cb.set_label('Median local correlation',fontsize=5.8,labelpad=1); cb.ax.tick_params(labelsize=5,length=2,pad=1); cb.outline.set_linewidth(.4)
panel(axC,'C',x=-.22,y=1.03)

axD=fig.add_subplot(gs[1,1]); d=D[D.record_type.eq('interfactor_lava')]
n=len(FACTOR_ORDER); mat=np.full((n,n),np.nan); prop=np.full((n,n),np.nan)
for _,r in d.iterrows():
    i=FACTOR_ORDER.index(r.factor1); j=FACTOR_ORDER.index(r.factor2)
    mat[i,j]=mat[j,i]=float(r.median_rho); prop[i,j]=prop[j,i]=float(r.proportion_fdr)
for i in range(n):
    for j in range(n):
        if i==j or np.isnan(mat[i,j]): continue
        axD.scatter(j,i,s=8+115*max(0,prop[i,j]),c=[mat[i,j]],cmap=CMAP,norm=Normalize(-1,1),edgecolor='none')
axD.set_xticks(range(n)); axD.set_xticklabels(FACTOR_ORDER,rotation=90,fontsize=5.3)
axD.set_yticks(range(n)); axD.set_yticklabels(FACTOR_ORDER,fontsize=5.3)
axD.set_xlim(-.6,n-.4); axD.set_ylim(n-.4,-.6); axD.tick_params(length=0)
axD.set_xticks(np.arange(-.5,n,1),minor=True); axD.set_yticks(np.arange(-.5,n,1),minor=True)
axD.grid(which='minor',color='#ECECEC',lw=.4)
for sp in axD.spines.values(): sp.set_visible(False)
for s,l in [(20,'0.1'),(65,'0.5'),(120,'1.0')]: axD.scatter([],[],s=s,color='#8A9BA8',label=l)
axD.legend(title='FDR-positive fraction',loc='upper left',bbox_to_anchor=(1.02,1),fontsize=5,title_fontsize=5.4,handletextpad=.2)
axD.text(1.02,.69,'UKB local-correlation\nresults only',transform=axD.transAxes,fontsize=5.1,color='#555555',va='top')
panel(axD,'D',x=-.12,y=1.03)

for ext in ['pdf','svg','png']:
    fig.savefig(ROOT/f'Figure4.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight',facecolor='white')
plt.close(fig)
