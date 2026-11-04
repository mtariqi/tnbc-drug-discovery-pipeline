from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).resolve().parent
x=pd.read_csv(ROOT/'expression_standardized.csv',index_col=0)
s=pd.read_csv(ROOT/'string_edges.tsv',sep='\t')
assert x.shape==(24,151),f'Unexpected matrix: {x.shape}'
assert not x.isna().any().any(), 'Missing protein abundances; use pairwise-complete n and review preprocessing'
string_pairs={tuple(sorted((str(a).strip().upper(),str(b).strip().upper()))) for a,b in zip(s.source,s.target)}
records=[]
for a,b in combinations(x.index,2):
    g1,g2=a.split('|')[-2],b.split('|')[-2]
    r,p=pearsonr(x.loc[a].to_numpy(dtype=float),x.loc[b].to_numpy(dtype=float))
    records.append({'gene1':g1,'gene2':g2,'pearson_r':r,'p_value':p,'n_samples':x.shape[1], 'string_edge':tuple(sorted((g1,g2))) in string_pairs})
df=pd.DataFrame(records)
df['bh_fdr']=multipletests(df.p_value,method='fdr_bh')[1]
df['significant']=df.bh_fdr<.05
df['abs_r']=df.pearson_r.abs()
df.sort_values('bh_fdr').to_csv(OUT/'all_276_correlations_verified.csv',index=False)
sig=df[df.significant].copy().sort_values('abs_r',ascending=False)
sig.to_csv(OUT/'significant_107_correlations_verified.csv',index=False)
matched=int(sig.string_edge.sum());unmatched=len(sig)-matched
assert len(df)==276 and len(sig)==107 and matched==44 and unmatched==63,(len(df),len(sig),matched,unmatched)
e=sig[((sig.gene1=='ERBB3')&(sig.gene2=='ABL1'))|((sig.gene1=='ABL1')&(sig.gene2=='ERBB3'))].iloc[0]
rank_all=int((sig.abs_r>abs(e.pearson_r)).sum()+1)
rank_unmatched=int((sig.loc[~sig.string_edge,'abs_r']>abs(e.pearson_r)).sum()+1)
print('tested',len(df),'significant',len(sig),'matched',matched,'unmatched',unmatched,'ERBB3-ABL1',e.to_dict(),'rank all',rank_all,'rank unmatched',rank_unmatched)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(13.2,8.4),facecolor='white')
gs=fig.add_gridspec(2,2,height_ratios=[1,1.18],width_ratios=[1,1.15],hspace=.44,wspace=.32,left=.09,right=.96,top=.88,bottom=.10)
navy='#17324D';blue='#3173A6';orange='#DA7B35';gray='#B9C4CC';green='#25847B'
fig.suptitle('CPTAC protein co-expression and STRING network concordance',x=.09,y=.975,ha='left',fontsize=16,fontweight='bold',color=navy)
fig.text(.09,.931,'Breast cancer proteomics · 24 kinases · 151 samples · Pearson correlations · BH-FDR < 0.05',fontsize=10,color='#536577')
ax=fig.add_subplot(gs[0,0]);ax.set_title('A  Statistical testing outcomes',loc='left',fontweight='bold',color=navy,pad=12)
parts=[44,63,169];colors=[blue,orange,gray];labels=['STRING-supported','No STRING edge','Not significant']
left=0
for n,c,l in zip(parts,colors,labels):
    ax.barh(0,n,left=left,color=c,height=.5,label=f'{l} (n={n})');ax.text(left+n/2,0,str(n),ha='center',va='center',color='white' if c!=gray else navy,fontsize=12,fontweight='bold');left+=n
ax.set_xlim(0,276);ax.set_ylim(-.6,.6);ax.set_yticks([]);ax.set_xlabel('Protein pairs tested (n = 276)');ax.legend(loc='lower left',bbox_to_anchor=(-.02,-.63),ncol=1,frameon=False,fontsize=9)
ax=fig.add_subplot(gs[0,1]);ax.set_title('B  STRING concordance among significant pairs',loc='left',fontweight='bold',color=navy,pad=12)
ax.bar(['STRING-supported','No STRING edge'],[44,63],color=[blue,orange],width=.52)
for i,n in enumerate([44,63]):ax.text(i,n+2,f'{n} ({n/107:.1%})',ha='center',fontsize=11,fontweight='bold',color=navy)
ax.set_ylim(0,78);ax.set_ylabel('BH-FDR significant pairs (n = 107)');ax.grid(axis='y',alpha=.14);ax.set_axisbelow(True)
ax=fig.add_subplot(gs[1,0]);ax.set_title('C  Distribution of observed protein correlations',loc='left',fontweight='bold',color=navy,pad=12)
for subset,c,l in [(df[~df.significant],gray,'Not significant'),(sig[sig.string_edge],blue,'STRING-supported'),(sig[~sig.string_edge],orange,'No STRING edge')]:
    ax.hist(subset.pearson_r,bins=np.linspace(-.8,.8,24),alpha=.75,color=c,label=f'{l} (n={len(subset)})',edgecolor='white',linewidth=.5)
ax.axvline(e.pearson_r,color=green,linestyle='--',linewidth=2,label='ERBB3–ABL1 (r = 0.506)')
ax.set_xlabel('Pearson correlation coefficient (r)');ax.set_ylabel('Pair count');ax.legend(frameon=False,fontsize=8,loc='upper left');ax.grid(axis='y',alpha=.12);ax.set_axisbelow(True)
ax=fig.add_subplot(gs[1,1]);ax.set_title('D  Highest |r| associations without a STRING edge',loc='left',fontweight='bold',color=navy,pad=12)
top=sig[~sig.string_edge].nlargest(8,'abs_r').copy().iloc[::-1];names=[f'{a}–{b}' for a,b in zip(top.gene1,top.gene2)]
barcols=[green if {a,b}=={'ERBB3','ABL1'} else orange for a,b in zip(top.gene1,top.gene2)]
ax.barh(names,top.pearson_r,color=barcols,height=.67)
for y,r in enumerate(top.pearson_r):ax.text(r+(.014 if r>=0 else -.014),y,f'{r:+.3f}',va='center',ha='left' if r>=0 else 'right',fontsize=8,color=navy)
ax.axvline(0,color='#536577',lw=.8);ax.set_xlim(-.73,.73);ax.set_xlabel('Pearson r');ax.grid(axis='x',alpha=.12);ax.set_axisbelow(True)
fig.text(.09,.025,f'ERBB3–ABL1: r = {e.pearson_r:.3f}; p = {e.p_value:.2e}; BH-FDR = {e.bh_fdr:.2e}; rank {rank_unmatched} among unmatched pairs by |r|.  Correlation does not establish direct binding.',fontsize=8.5,color='#536577')
for ext in ('png','svg','pdf'):fig.savefig(OUT/f'Figure5_CPTAC_STRING_VERIFIED.{ext}',dpi=400,facecolor='white')
plt.close(fig)
summary=pd.DataFrame([{'samples':151,'kinases':24,'tested_pairs':276,'significant_bh_fdr_005':107,'string_supported':44,'not_in_string':63,'erbb3_abl1_r':e.pearson_r,'erbb3_abl1_p':e.p_value,'erbb3_abl1_bh_fdr':e.bh_fdr,'erbb3_abl1_rank_among_unmatched_absr':rank_unmatched}]);summary.to_csv(OUT/'summary_verified.csv',index=False)
