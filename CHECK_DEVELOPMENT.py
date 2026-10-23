"""Verify archived weight screening and shared endpoint identities."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent;D=ROOT/'development';R=D/'records'
for row in pd.read_csv(D/'source_index.csv').itertuples():
    h=hashlib.sha256((D/row.file).read_bytes()).hexdigest()
    assert h==row.sha256==row.source_sha256,row.file
for row in pd.read_csv(D/'stage_index.csv').itertuples():
    assert all((D/name).is_file() for name in row.records.split('|'))
data=pd.read_csv(R/'lambda_screening_endpoints225.csv')
assert len(data)==225 and not data.duplicated(['benchmark','config','repeat']).any()
assert set(data.repeat)=={101,102,103}
means=data.groupby(['config','benchmark']).common_rmse.mean().unstack()
ratios=means/means.loc['l25c20']
ranking=pd.DataFrame({'geometric_mean_ratio':np.exp(np.log(ratios).mean(axis=1)),'worst_ratio':ratios.max(axis=1),'GL_ratio':ratios.GL1D,'Peaks_ratio':ratios.Peaks2D,'H3_ratio':ratios.H3_historical})
stored=pd.read_csv(R/'lambda_screening_ranking25.csv').set_index('config')
np.testing.assert_allclose(ranking.loc[stored.index,stored.columns],stored,rtol=1e-12,atol=1e-14)
ranked=ranking.reset_index().sort_values(['geometric_mean_ratio','worst_ratio','config'])
selection=json.loads((R/'lambda_selection.json').read_text())
assert ranked.iloc[0]['config']==selection['candidate']=='l35c10'
base=pd.read_csv(ROOT/'data/numerical/baseline_endpoints540.csv')
base=base[base.version=='FDCVV'].set_index(['benchmark','repeat']).sort_index()
assert len(base)==90
checks=[]
for file,version in [('K_schedule_pilot_endpoints210.csv','K432_33_50'),('K_schedule_extend15_candidate45.csv','K432_33_50'),('K_schedule_extend30_endpoints360.csv','Full')]:
    frame=pd.read_csv(R/file);g=frame[frame.version==version].set_index(['benchmark','repeat']).sort_index()
    assert not g.index.duplicated().any()
    np.testing.assert_allclose(g[['common_rmse','common_r2']],base.loc[g.index,['common_rmse','common_r2']],rtol=1e-12,atol=1e-14)
    checks.append(dict(file=file,matched_endpoints=len(g)))
reuse=pd.read_csv(D/'trajectory_reuse90.csv').set_index(['benchmark','repeat']).sort_index()
for file,field,label in [('endpoints720.csv','version','Full'),('lambda_endpoints2250.csv','config',reuse.grid_config.iloc[0])]:
    a=pd.read_csv(ROOT/'data/numerical'/file)
    a=a[a[field]==label].set_index(['benchmark','repeat']).sort_index()
    assert a.index.equals(base.index)
    np.testing.assert_allclose(a[['common_rmse','common_r2','true_nrmse']],base[['common_rmse','common_r2','true_nrmse']],rtol=1e-12,atol=1e-14)
followup=pd.read_csv(R/'lambda_followup_endpoints45.csv')
assert len(followup)==45 and set(followup.repeat)==set(range(106,111))
report=dict(status='PASS',weight_screening_rows=len(data),weight_candidate=selection['candidate'],score=float(ranked.iloc[0].geometric_mean_ratio),followup_rows=len(followup),shared_final_endpoints=90,developmental_subsets=checks,new_model_fits=0)
print(json.dumps(report,indent=2))
