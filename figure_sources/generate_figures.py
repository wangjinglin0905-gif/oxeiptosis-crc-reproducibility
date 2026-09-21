from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from PIL import Image

O=Path(__file__).resolve().parents[1]
T=O/'supporting_information';F=O/'figures';V=O/'figure_sources'
F.mkdir(exist_ok=True);V.mkdir(exist_ok=True)
def rows(n):return list(csv.DictReader((T/n).open(encoding='utf-8-sig')))
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':9,'legend.fontsize':8.5,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'savefig.facecolor':'white'})
blue='#2166A5';orange='#C26A13';gray='#707070';coh=['TCGA-COADREAD','GSE39582','GSE17538'];manifest=[]
def save(fig,k):
    fig.canvas.draw();ren=fig.canvas.get_renderer();outside=[]
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text():
            b=t.get_window_extent(ren)
            if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.x1+.5 or b.y1>fig.bbox.y1+.5:outside.append(t.get_text())
    for ext in ['svg','pdf','png']:fig.savefig(V/f'Fig{k}.{ext}',dpi=300)
    p=F/f'Fig{k}.tif';fig.savefig(p,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    with Image.open(p) as im:rgb=im.convert('RGB')
    rgb.save(p,dpi=(600,600),compression='tiff_lzw')
    manifest.append({'figure':k,'width_inches':fig.get_figwidth(),'height_inches':fig.get_figheight(),'pixels':rgb.size,'dpi':600,'mode':'RGB','outside_text':outside,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    plt.close(fig)
r=[x for x in rows('S7_OS_all45_R.csv') if x['score']=='OX4']
fig=plt.figure(figsize=(6.5,4.25));ax=fig.add_axes([.30,.15,.30,.79]);tx=fig.add_axes([.63,.15,.36,.79]);tx.axis('off');tx.set_ylim(-.6,10.7)
ys=[9.5,8.6,7.7,6,5.1,4.2,2.5,1.6,.7];i=0
for co in coh:
    z=[x for x in r if x['cohort']==co]
    ax.text(-.83,ys[i]+.62,f"{co}  {z[0]['n']}/{z[0]['events']}",transform=ax.get_yaxis_transform(),fontsize=9,fontweight='bold')
    for mo in ['M0','M1','M2']:
        x=next(x for x in z if x['model']==mo);hr,lo,hi=map(float,[x['HR'],x['lo95'],x['hi95']]);y=ys[i];i+=1;c=blue if mo=='M1' else gray
        ax.errorbar(hr,y,xerr=[[hr-lo],[hi-hr]],fmt='s' if mo=='M1' else 'o',color=c,markersize=4,capsize=2,lw=1.2)
        ax.text(-.1,y,mo,transform=ax.get_yaxis_transform(),ha='right',va='center',color=c)
        tx.text(0,y,f'{hr:.3f} ({lo:.3f}–{hi:.3f})',va='center',fontsize=9,color=c)
ax.set_ylim(-.6,10.7);ax.set_xscale('log');ax.set_xlim(.47,1.62);ax.set_xticks([.5,.75,1,1.5],labels=['0.5','0.75','1','1.5']);ax.xaxis.set_minor_locator(NullLocator());ax.set_yticks([]);ax.spines['left'].set_visible(False);ax.axvline(1,color=gray,ls='--',lw=.8);ax.set_xlabel('Hazard ratio per reference SD');fig.text(.64,.96,'HR (95% CI)',fontsize=9,fontweight='bold');save(fig,1)
b=rows('S24_matched_null_all6000_R.csv');sm=rows('S23_matched_null_summary.csv');mx=max(abs(np.log(float(x['HR']))) for x in b)
fig,axs=plt.subplots(1,3,figsize=(6.5,3.0),sharey=True);fig.subplots_adjust(left=.10,right=.98,bottom=.26,top=.83,wspace=.22)
for i,(ax,co) in enumerate(zip(axs,coh)):
    for de,c,label in [('M1_msd',blue,'MSD'),('M2_msd_coexpr',orange,'MSD + OX4 correlation')]:
        v=np.sort([abs(np.log(float(x['HR']))) for x in b if x['cohort']==co and x['design']==de]);ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',color=c,label=label,lw=1.3)
    ob=float(next(x['observed_abs_logHR'] for x in sm if x['cohort']==co));ax.axvline(ob,color='#222222',ls='--',lw=1);ax.set_xlim(0,mx*1.03);ax.set_ylim(0,1.03);ax.set_xticks([0,.2,.4]);ax.set_xlabel('Absolute log HR');ax.set_title(co);ax.text(-.12,1.12,chr(97+i),transform=ax.transAxes,fontweight='bold',fontsize=11)
axs[0].set_ylabel('Cumulative fraction');fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.54,0),ncol=2,frameon=False);save(fig,2)
c=rows('S9_transcript_and_subcohort_correlations.csv');refs=['S_mito_c','S_prolif_c','S_nrf2'];labs=['MitoCarta','Proliferation','NRF2 targets']
fig,axs=plt.subplots(1,2,figsize=(7.3,3.2),gridspec_kw={'width_ratios':[1,1.05]});fig.subplots_adjust(left=.12,right=.98,bottom=.33,top=.84,wspace=.58)
v=np.array([[float(next(x['rho'] for x in c if x['cohort']==co and x['reference']==rf)) for co in ['TCGA','GSE39582','GSE17538']] for rf in refs]);axs[0].imshow(v,vmin=0,vmax=1,cmap='Blues',aspect='auto')
for i in range(3):
    for j in range(3):axs[0].text(j,i,f'{v[i,j]:.3f}',ha='center',va='center',color='white' if v[i,j]>.65 else '#222222')
axs[0].set_yticks(range(3),labs);axs[0].set_xticks(range(3),coh,rotation=25,ha='right');axs[0].tick_params(length=0)
ax=axs[1]
for i,rf in enumerate(refs):
    x=next(x for x in c if x['cohort']=='GSE17538' and x['reference']==rf);raw,adj=float(x['rho']),float(x['subcohort_adjusted_rho']);ax.plot([adj,raw],[i,i],color=gray,lw=1);ax.plot(raw,i,'o',color=blue,label='Raw' if i==0 else None);ax.plot(adj,i,'D',color=orange,label='Subcohort adjusted' if i==0 else None)
ax.set_yticks(range(3),labs);ax.invert_yaxis();ax.set_xlim(.4,.86);ax.set_xticks([.4,.6,.8]);ax.set_xlabel('Rank correlation');ax.set_title('GSE17538');fig.legend(*ax.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.55,-.01),ncol=2,frameon=False)
for i,ax in enumerate(axs):ax.text(-.15,1.12,chr(97+i),transform=ax.transAxes,fontweight='bold',fontsize=11)
save(fig,3)
p=rows('S11_reference_pair_all108.csv');ds=rows('S12_reference_pairs_by_cohort.csv')
fig,axs=plt.subplots(1,3,figsize=(6.5,3.0),sharey=True);fig.subplots_adjust(left=.11,right=.98,bottom=.20,top=.79,wspace=.22)
for i,(ax,co) in enumerate(zip(axs,coh)):
    g=[x for x in p if x['cohort']==co];ax.scatter([float(x['jaccard_measured_clean_audit']) for x in g],[float(x['rho_clean_audit']) for x in g],s=18,color=blue,alpha=.7,edgecolors='white',linewidth=.3)
    ax.axvline(.05,color=gray,ls='--',lw=.7);ax.axhline(.5,color=gray,ls='--',lw=.7);ax.set_xlim(-.01,.2);ax.set_ylim(-.2,1.05);ax.set_xticks([0,.1,.2],labels=['0','0.10','0.20']);ax.set_xlabel('Jaccard overlap');n=next(x['n'] for x in ds if x['cohort']==co);ax.set_title(f'{co}\nn = {n}');ax.text(-.12,1.20,chr(97+i),transform=ax.transAxes,fontweight='bold',fontsize=11)
axs[0].set_ylabel('Reference-score Spearman rho');save(fig,4)
(O/'figure_sources/figure_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps(manifest,indent=2))
