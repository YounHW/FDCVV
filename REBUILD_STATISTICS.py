"""Reconstruct numerical summaries from saved run-level RMSE; no model fitting."""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
import pandas as pd
import scipy
from scipy.io import loadmat
from scipy.stats import ttest_rel

ROOT=Path(__file__).resolve().parent
BEN=['GL1D','Peaks2D','H3_historical']
CONTROLS=['RandomRegions','NoVoronoi','FixedK2','FixedK4','FixedLambda','VarianceOnly','DistanceOnly']
METRICS=['common_rmse','true_nrmse','common_r2']

def read(name):
    return pd.read_csv(ROOT/'data/numerical'/name)

def summarize(df,keys,metrics=METRICS):
    rows=[]
    for labels,g in df.groupby(keys,sort=True):
        row=dict(zip(keys,labels if isinstance(labels,tuple) else (labels,)))
        row['n']=len(g)
        for metric in metrics:
            a=g[metric].to_numpy()
            row.update({metric+'_'+k:float(v) for k,v in dict(mean=a.mean(),std=a.std(ddof=1),median=np.median(a),min=a.min(),max=a.max()).items()})
            row[metric+'_sem']=float(a.std(ddof=1)/np.sqrt(len(a)))
        rows.append(row)
    return pd.DataFrame(rows)

def paired_effects(end):
    rng=np.random.default_rng(20261013);rows=[]
    for b in BEN:
        p=end[end.benchmark==b].pivot(index='repeat',columns='version',values='true_nrmse').sort_index()
        assert list(p.index)==list(range(106,136))
        full=p.Full.to_numpy();indices=rng.integers(0,30,size=(20000,30))
        for version in CONTROLS:
            control=p[version].to_numpy();delta=control-full
            boot_rel=100*(control[indices].mean(axis=1)/full[indices].mean(axis=1)-1)
            boot_delta=delta[indices].mean(axis=1);test=ttest_rel(control,full)
            signs=rng.choice([-1.,1.],size=(100000,30))
            ps=(1+np.count_nonzero(np.abs((signs*delta).mean(axis=1))>=abs(delta.mean())))/100001
            lo,hi=np.quantile(boot_rel,[.025,.975]);dlo,dhi=np.quantile(boot_delta,[.025,.975])
            rows.append(dict(benchmark=b,version=version,n=30,full_mean=full.mean(),control_mean=control.mean(),difference=delta.mean(),difference_ci_low=dlo,difference_ci_high=dhi,relative_pct=100*(control.mean()/full.mean()-1),relative_ci_low=lo,relative_ci_high=hi,t=float(test.statistic),df=29,p=float(test.pvalue),p_signflip=ps,full_wins=int((delta>1e-10).sum()),ties=int((abs(delta)<=1e-10).sum()),full_losses=int((delta < -1e-10).sum())))
    result=pd.DataFrame(rows)
    for col in ['p','p_signflip']:
        order=np.argsort(result[col].to_numpy())
        result.loc[result.index[order],col+'_holm21']=np.minimum(1,np.maximum.accumulate(result[col].to_numpy()[order]*(21-np.arange(21))))
    result['metric']='true_nrmse'
    return result

def grid_gaps(summary):
    rows=[]
    for metric in ['common_rmse','true_nrmse']:
        field=metric+'_mean';minimum=summary.groupby('benchmark')[field].min()
        for config,g in summary.groupby('config',sort=True):
            row=dict(metric=metric,config=config,lambda0=g.lambda0.iloc[0],p_cross=g.p_cross.iloc[0])
            delta={b:100*(g.loc[g.benchmark==b,field].iloc[0]/minimum[b]-1) for b in BEN}
            row.update({b+'_gap_pct':v for b,v in delta.items()})
            row.update(worst_gap_pct=max(delta.values()),mean_gap_pct=np.mean(list(delta.values())))
            rows.append(row)
    result=pd.DataFrame(rows)
    result['descriptive_rank']=result.groupby('metric').worst_gap_pct.rank(method='min')
    return result

def stopping_summary(outcomes):
    rng=np.random.default_rng(26101405);rows=[]
    for b in BEN:
        for threshold in [.005,.01,.02]:
            g=outcomes[(outcomes.benchmark==b)&(outcomes.threshold==threshold)].sort_values('repeat')
            assert list(g.repeat)==list(range(106,136))
            v=lambda key:g[key].to_numpy(dtype=float)
            indices=rng.integers(0,30,(20000,30));saved=v('saved_points')
            lo,hi=np.quantile(saved[indices].mean(axis=1),[.025,.975])
            row=dict(benchmark=b,threshold=threshold,n=30,early_count=int(sum(saved>0)),triggered_count=int(sum(v('triggered'))),mean_saved=saved.mean(),saved_ci_low=lo,saved_ci_high=hi,mean_stop_step=v('stop_step').mean(),mean_training_n=v('training_n').mean())
            for metric in ['rmse','common_nrmse','r2','true_nrmse']:
                a=v('stop_'+metric);full=v('full_'+metric);delta=a-full
                lo,hi=np.quantile(delta[indices].mean(axis=1),[.025,.975])
                row.update({metric+'_stop_mean':a.mean(),metric+'_full_mean':full.mean(),metric+'_difference':delta.mean(),metric+'_difference_ci_low':lo,metric+'_difference_ci_high':hi})
                if metric!='r2':
                    rel=100*(a[indices].mean(axis=1)/full[indices].mean(axis=1)-1)
                    lo,hi=np.quantile(rel,[.025,.975])
                    row.update({metric+'_relative_pct':100*(a.mean()/full.mean()-1),metric+'_relative_ci_low':lo,metric+'_relative_ci_high':hi})
            rows.append(row)
    return pd.DataFrame(rows)

def check_table(actual,stored,keys,columns=None):
    assert not actual.duplicated(keys).any() and not stored.duplicated(keys).any()
    a=actual.sort_values(keys).reset_index(drop=True);b=stored.sort_values(keys).reset_index(drop=True)
    assert len(a)==len(b)
    for col in columns or list(stored.columns):
        if pd.api.types.is_numeric_dtype(b[col]):
            np.testing.assert_allclose(a[col],b[col],rtol=1e-9,atol=1e-11,equal_nan=True,err_msg=col)
        else:
            assert a[col].fillna('').equals(b[col].fillna('')),col

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/rebuilt_statistics')
    parser.add_argument('--check',action='store_true',help='Compare reconstructed columns/statistics with the supplied tables')
    args=parser.parse_args();out=args.output.resolve()
    for folder in ['data','tests','development']:
        assert not out.is_relative_to((ROOT/folder).resolve()),'Use a separate output directory'
    out.mkdir(parents=True,exist_ok=True)
    ranges={};bindings={}
    for b in BEN:
        p=ROOT/'tests'/f'{b}_test.mat';truth=np.asarray(loadmat(p)['truth']).ravel()
        assert np.isfinite(truth).all() and np.ptp(truth)>0
        ranges[b]=float(np.ptp(truth));bindings[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    tables={};checked=[]
    for name,key,count in [('baseline_endpoints540.csv',['benchmark','version','repeat'],540),('baseline_curves16200.csv',['benchmark','version','repeat','step'],16200),('endpoints720.csv',['benchmark','version','repeat'],720),('lambda_endpoints2250.csv',['benchmark','config','repeat'],2250)]:
        original=read(name);df=original.copy();assert len(df)==count and not df.duplicated(key).any()
        df['test_response_range']=df.benchmark.map(ranges);df['true_nrmse']=df.common_rmse/df.test_response_range
        assert np.isfinite(df.true_nrmse).all()
        tables[name]=df
        if args.check:
            check_table(df,original,key,columns=key+['test_response_range','true_nrmse']);checked.append(name)
    baseline=tables['baseline_endpoints540.csv'];curves=tables['baseline_curves16200.csv']
    end=tables['endpoints720.csv'];grid=tables['lambda_endpoints2250.csv']
    check_table(curves[curves.step==30],baseline,['benchmark','version','repeat'],columns=['benchmark','version','repeat']+METRICS)
    full=baseline[baseline.version=='FDCVV'].drop(columns='version')
    check_table(end[end.version=='Full'],full,['benchmark','repeat'],columns=['benchmark','repeat']+METRICS)
    retained=grid[np.isclose(grid.lambda0,.35)&np.isclose(grid.p_cross,.1)]
    check_table(retained,full,['benchmark','repeat'],columns=['benchmark','repeat']+METRICS)
    tables['baseline_summary18.csv']=summarize(baseline,['benchmark','version'])
    tables['baseline_curve_summary540.csv']=summarize(curves,['benchmark','version','step','samples'])
    tables['ablation_summary24.csv']=summarize(end,['benchmark','version'])
    tables['paired_effects21.csv']=paired_effects(end)
    ls=summarize(grid,['config','benchmark','lambda0','p_cross'],['common_rmse','common_nrmse','common_r2','true_nrmse'])
    tables['lambda_summary75.csv']=ls
    tables['lambda_cross_benchmark_gap50.csv']=grid_gaps(ls)
    stops=read('stopping_current_outcomes270.csv').copy()
    for prefix in ['stop','full']:
        stops[prefix+'_true_nrmse']=stops[prefix+'_rmse']/stops.benchmark.map(ranges)
    tables['stopping_current_outcomes270.csv']=stops
    tables['stopping_summary9.csv']=stopping_summary(stops)
    if args.check:
        for name,keys in [('paired_effects21.csv',['benchmark','version']),('lambda_summary75.csv',['config','benchmark']),('lambda_cross_benchmark_gap50.csv',['metric','config']),('stopping_current_outcomes270.csv',['benchmark','repeat','threshold']),('stopping_summary9.csv',['benchmark','threshold'])]:
            check_table(tables[name],read(name),keys);checked.append(name)
    for name,df in tables.items():
        df.to_csv(out/name,index=False)
    for name in ['baseline_endpoints540.csv','baseline_curves16200.csv','endpoints720.csv','lambda_endpoints2250.csv','stopping_current_outcomes270.csv']:
        p=ROOT/'data/numerical'/name;bindings[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    effects=tables['paired_effects21.csv'];sig=effects.p_holm21<.05
    report=dict(status='PASS',saved_table_checks=checked,test_response_ranges=ranges,input_sha256=bindings,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,bootstrap=dict(ablation_seed=20261013,stopping_seed=26101405,draws=20000,quantile_method='linear',interval='pointwise paired percentile 95%'),signflip=dict(seed=20261013,draws=100000,pseudocount=1),holm_family_size=21,classification=dict(full_better=int((sig&(effects.difference>0)).sum()),control_better=int((sig&(effects.difference<0)).sum()),not_significant=int((~sig).sum())),new_model_fits=0)
    (out/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='PASS',checked_tables=len(checked),exported_tables=len(tables),classification=report['classification'],output=str(out)),indent=2))

if __name__=='__main__':
    main()
