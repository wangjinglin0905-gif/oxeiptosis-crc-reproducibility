"""Portable audit of supplied morphology tables, not a raw-count pipeline.
Usage: python verify_frozen_morphology.py
Dependencies: numpy, pandas. Writes nothing; exits nonzero on discrepancy.
"""
from pathlib import Path
import itertools,json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1]/'tables'
def read(name):
    paths=list(R.rglob(name));assert len(paths)==1,(name,len(paths))
    return pd.read_csv(paths[0])
sec=read('S33_section_results.csv');p=read('S32_patient_effects.csv');inf=read('S32_patient_inference.csv');b=read('S36_all_construct_summaries.csv');pos=read('S34_matched_positions.csv');spl=read('S37_technical_splits.csv')
def equal(a,b):assert np.allclose(a,b,rtol=0,atol=1e-10,equal_nan=True),(a,b)
cols=['cliff_delta','mean_difference_sd','depth_adjusted_delta','median_split_rho']
calc=sec.groupby(['patient','construct'])[cols].mean().sort_index()
saved=p.set_index(['patient','construct']).sort_index()[cols]
equal(calc.to_numpy(),saved.to_numpy())
signs=np.array(list(itertools.product([-1,1],repeat=6)))
ps=[]
for k,c in enumerate(['OX4','epithelial','stromal']):
    x=p[p.construct==c].sort_values('patient').cliff_delta.to_numpy();assert len(x)==6
    r=inf[inf.construct==c].iloc[0]
    z=x[np.random.default_rng(20261118+k).integers(0,6,(10000,6))].mean(axis=1)
    pv=np.mean(abs((signs*x).mean(axis=1))>=abs(x.mean())-1e-12);ps.append(pv)
    equal([x.mean(),*np.quantile(z,[.025,.975]),pv],[r.mean_delta,r.bootstrap_lo,r.bootstrap_hi,r.exact_signflip_p])
ps=np.array(ps);order=np.argsort(ps);qs=np.empty(3);qs[order]=np.minimum(1,np.minimum.accumulate((ps[order]*3/np.arange(1,4))[::-1])[::-1])
equal(qs,inf.set_index('construct').loc[['OX4','epithelial','stromal'],'BH_q_three_scores'].to_numpy())
assert len(sec)==36 and len(spl)==37080 and len(b)==1236
for _,r in pos.iterrows():
    d=b[b['sample']==r['sample']];ox=d[d.construct=='OX4'].iloc[0];null=d[d.construct.str.startswith('random_')];assert len(null)==100
    equal([100*(abs(null.cliff_delta)<=abs(ox.cliff_delta)).mean(),100*(null.median_split_rho<=ox.median_split_rho).mean()],[r.abs_delta_percentile,r.split_rho_percentile])
for (sid,c),d in spl.groupby(['sample','construct']):
    assert len(d)==30
    equal(d.rho_all.median(),b[(b['sample']==sid)&(b.construct==c)].iloc[0].median_split_rho)
print(json.dumps({'status':'passed','patients':6,'sections':12,'construct_summaries':len(b),'technical_records':len(spl),'numpy':np.__version__,'pandas':pd.__version__,'scope':'frozen-table aggregation, bootstrap, sign-flip, BH and matched-position verification; original scoring/count splitting not rerun'}))
