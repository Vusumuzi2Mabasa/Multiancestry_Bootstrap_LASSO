"""Audit complete paired results and summarize fixed-stratum Monte Carlo error."""
from common import *
import pandas as pd,math
from statistics import NormalDist
METRICS=['any_selection','fdp','recall','ld_recall_ceu','ld_recall_african','mse_ceu','mse_african','specific_estimate','seconds','algorithmic_jaccard','perturb_jaccard','B200_jaccard']
def aggregate(g,col):
    cells=[x[col].dropna().astype(float).values for _,x in g.groupby(['region','population'])]
    cells=[x for x in cells if len(x)]
    if not cells:return (float('nan'),float('nan'),0)
    mean=np.mean([x.mean() for x in cells]);se=np.sqrt(sum(np.var(x,ddof=1)/len(x) if len(x)>1 else 0 for x in cells))/len(cells)
    return float(mean),float(se),sum(map(len,cells))
def main():
    files=list((ROOT/'results/validation').glob('*/metrics.tsv'));assert len(files)==12
    d=pd.concat([pd.read_csv(p,sep='\t',keep_default_na=False,na_values=['']) for p in files],ignore_index=True)
    variants=pd.concat([pd.read_csv(p,sep='\t',keep_default_na=False,na_values=['']) for p in (ROOT/'results/validation').glob('*/variants.tsv')],ignore_index=True)
    key=['region','population','scenario','replicate','ld','method']
    assert not d.duplicated(key).any();assert len(d)==12*180*2*4
    assert not variants.duplicated(key+['index']).any() and len(variants)==len(d)*8
    assert np.isfinite(d[['mse_ceu','mse_african','seconds']]).all().all()
    assert d.groupby(key[:-1]).size().eq(4).all()
    # Reconstruct metrics independently from the per-variant table.
    max_mse_error=0
    index=d.set_index(key)
    cutoff_cache={p.parent.name:json.loads(p.read_text())['cutoffs'] for p in (ROOT/'results/validation').glob('*/thresholds.json')}
    for values,g in variants.groupby(key):
        m=index.loc[values];truth=(g.truth_ceu!=0)|(g.truth_african!=0);sel=g.selected==1
        assert sel.sum()==m.selected and (sel&truth).sum()==m.tp and (sel&~truth).sum()==m.fp
        assert abs(m.fdp-(sel&~truth).sum()/max(1,sel.sum()))<1e-12
        for pop in ['ceu','african']:
            mse=float(((g['estimate_'+pop]-g['truth_'+pop])**2).mean());max_mse_error=max(max_mse_error,abs(m['mse_'+pop]-mse))
        if values[-1]=='Bootstrap_LASSO_Manc':
            cut=cutoff_cache[f'r{values[0]}_{values[1]}'][values[-2]][values[-1]]
            assert int((g.frequency>cut).sum())==m.proposed
    assert max_mse_error<1e-12
    freeze=json.loads((ROOT/'protocol/validation_threshold_freeze.json').read_text());assert freeze['code_sha256']==codehash()
    for p,h in freeze['thresholds'].items():assert sha(ROOT/p)==h
    complete=[json.loads(p.read_text()) for p in (ROOT/'results/validation').glob('*/*_complete.json')];assert len(complete)==24
    assert all(x['code_sha256']==codehash() and x['failed_fits']==0 for x in complete)
    summary=[]
    for (ld,scenario,method),g in d.groupby(['ld','scenario','method']):
        row=dict(ld=ld,scenario=scenario,method=method,n=len(g))
        for metric in METRICS:
            mean,se,n=aggregate(g,metric);row[metric]=mean;row[metric+'_mcse']=se
        row['rmse_african']=math.sqrt(row['mse_african']);row['rmse_ceu']=math.sqrt(row['mse_ceu'])
        if scenario=='null':
            rate=row['any_selection'];n=len(g);z=1.96;center=(rate+z*z/(2*n))/(1+z*z/n);half=z*math.sqrt(rate*(1-rate)/n+z*z/(4*n*n))/(1+z*z/n)
            row['null_wilson_lower']=center-half;row['null_wilson_upper']=center+half
        else:row['null_wilson_lower']=row['null_wilson_upper']=float('nan')
        summary.append(row)
    table(ROOT/'results/validation_summary.tsv',summary)
    paired=[]
    primary=d[d.method=='Bootstrap_LASSO_Manc']
    for method in ['Manc_stepwise','Exact_rho08','MultiSuSiE_L5']:
        comp=d[d.method==method];join=primary.merge(comp,on=key[:-1],suffixes=('_a','_b'),validate='one_to_one')
        for (ld,sc),g in join.groupby(['ld','scenario']):
            for metric in METRICS:
                vals=g[['region','population']].copy();vals['difference']=g[metric+'_a']-g[metric+'_b']
                est,se,n=aggregate(vals,'difference')
                if n:paired.append(dict(ld=ld,scenario=sc,comparator=method,metric=metric,n=n,difference=est,mcse=se,lower=est-1.96*se,upper=est+1.96*se))
    table(ROOT/'results/paired_differences.tsv',paired)
    cuts=[]
    for p in (ROOT/'results/validation').glob('*/thresholds.json'):
        v=json.loads(p.read_text())
        for ld,methods in v['cutoffs'].items():
            for method,c in methods.items():cuts.append(dict(stratum=p.parent.name,ld=ld,method=method,cutoff=c))
    table(ROOT/'results/threshold_summary.tsv',cuts)
    table(ROOT/'results/fit_status.tsv',d.groupby(['ld','method','status']).size().reset_index(name='count').to_dict('records'))
    save(ROOT/'results/audit.json',dict(passed=True,unique_evaluation_datasets=12*180,unique_calibration_datasets=12*199,paired_method_rows=len(d),variant_rows=len(variants),failed_fits=0,nonconverged=int((d.converged==0).sum()),
      max_kkt=max(x['max_kkt'] for x in complete),max_mse_reconstruction_error=max_mse_error,code_sha256=codehash(),
      total_joint_drops=int(d.joint_dropped.sum()),causal_filtered=int((d.causal_count-d.retained_causal).sum()),
      result_hashes={str(p.relative_to(ROOT)):sha(p) for p in files}))
    plots(pd.DataFrame(summary),pd.DataFrame(paired),d)
    print('Independent reconstruction audit passed;',len(d),'paired method rows')
def plots(s,p,d):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':14,'axes.titlesize':14,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':180})
    names={'Bootstrap_LASSO_Manc':'Bootstrap LASSO + Manc joint','Manc_stepwise':'Original Manc stepwise','Exact_rho08':'Exact Bayes rho 0.8','MultiSuSiE_L5':'MultiSuSiE L5'}
    colors={'Bootstrap_LASSO_Manc':'#007F86','Manc_stepwise':'#BD5436','Exact_rho08':'#696969','MultiSuSiE_L5':'#8F69A8'}
    primary=['Bootstrap_LASSO_Manc','Manc_stepwise'];scenarios=['one','shared','heterogeneous','three']
    fig,axs=plt.subplots(1,3,figsize=(12,4.8));a=axs[0]
    for i,m in enumerate(primary):
        for j,ld in enumerate(['matched','in_sample']):
            r=s[(s.method==m)&(s.scenario=='null')&(s.ld==ld)].iloc[0]
            a.errorbar(i+(j-.5)*.18,r.any_selection,yerr=1.96*r.any_selection_mcse,fmt='o' if j==0 else 's',color=colors[m],alpha=1 if j==0 else .45,capsize=3)
    a.axhline(.05,color='black',ls=':',lw=1);a.set_xticks([0,1],['Bootstrap','Stepwise']);a.set_title('Held-out global-null error');a.set_ylabel('Probability of any selection');a.set_ylim(bottom=0)
    for a,metric,title in zip(axs[1:],['recall','fdp'],['Exact recovery','False discovery proportion']):
        for m in primary:
            for ld in ['matched','in_sample']:
                sub=s[(s.method==m)&(s.ld==ld)].set_index('scenario').loc[scenarios]
                a.errorbar(range(4),sub[metric],yerr=1.96*sub[metric+'_mcse'],color=colors[m],ls='-' if ld=='matched' else '--',alpha=1 if ld=='matched' else .45,marker='o',capsize=2,label=names[m] if ld=='matched' else None)
        a.set_xticks(range(4),['One','Shared','Het.','Three']);a.set_title(title);a.set_ylim(bottom=0)
    handles,labels=axs[1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.01),frameon=False,ncol=2,fontsize=13)
    fig.suptitle('SIMULATED phenotypes | solid: external LD; dashed/squares: in-sample LD',fontsize=13)
    fig.subplots_adjust(left=.08,right=.99,bottom=.25,top=.8,wspace=.48);fig.savefig(ROOT/'figures/primary_comparison.png',bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(12,4.8))
    for a,metric,title in zip(axs,['recall','fdp','rmse_african'],['Exact recall','Mean FDP','African-source RMSE']):
        for m in names:
            sub=s[(s.method==m)&(s.ld=='matched')].set_index('scenario').loc[scenarios]
            a.plot(range(4),sub[metric],marker='o',color=colors[m],label=names[m]);a.set_xticks(range(4),['One','Shared','Het.','Three'])
        a.set_title(title);a.set_ylim(bottom=0)
    handles,labels=axs[1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.005),frameon=False,ncol=2,fontsize=13)
    fig.suptitle('SIMULATED phenotypes | secondary Bayesian comparators use posterior means',fontsize=13)
    fig.subplots_adjust(left=.08,right=.99,bottom=.32,top=.8,wspace=.48);fig.savefig(ROOT/'figures/secondary_comparison.png',bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,4.5))
    for m in primary:
        sub=s[(s.method==m)&(s.ld=='matched')].set_index('scenario').loc[['null']+scenarios]
        axs[0].plot(range(5),sub.perturb_jaccard,'o-',color=colors[m],label=names[m]);axs[1].plot(range(5),sub.seconds,'o-',color=colors[m])
    for a in axs:a.set_xticks(range(5),['Null','One','Shared','Het.','Three'])
    axs[0].set_title('Common summary perturbation stability');axs[0].set_ylabel('Mean Jaccard');axs[0].set_ylim(0,1.04)
    axs[1].set_title('Measured pipeline runtime');axs[1].set_ylabel('Seconds per context')
    axs[0].legend(frameon=False,fontsize=11)
    fig.suptitle('SIMULATED phenotypes | matched LD | empty-empty Jaccard = 1',fontsize=13)
    fig.tight_layout();fig.savefig(ROOT/'figures/stability_runtime.png',bbox_inches='tight');plt.close(fig)
if __name__=='__main__':main()
