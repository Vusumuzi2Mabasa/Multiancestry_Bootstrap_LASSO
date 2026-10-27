"""Paired simulation; calibration is a separate completed stage before evaluation."""
from common import *
from lasso_bridge import run_batch,case,frequencies
from manc import run as manc_run,entry_score
from bayes import custom,multisusie,match_recovery
import argparse,time,math,copy
from concurrent.futures import ProcessPoolExecutor,as_completed
SCENARIOS=['null','shared','heterogeneous','one','three']
def seed(phase,region,pop,scenario,rep,offset=0):
    sid=['calibration']+SCENARIOS
    return CFG['bootstrap_seed_offsets'][phase]+region*100000+CFG['populations'].index(pop)*10000+sid.index(scenario)*1000+rep+offset
def lambdas(phase):
    return CFG['lambda_grid'] if phase=='tuning' else [json.loads((ROOT/'protocol/tuning_freeze.json').read_text())['lambda']]
def label(lam,phase):return f'Bootstrap_lambda{lam:g}' if phase=='tuning' else 'Bootstrap_LASSO_Manc'
def bayesfit(stats,rs,keep,method):
    # Same eight-variant harmonized inputs. The validated exact comparator requires all 8.
    # Ineligible markers have their selection scores masked; latent model includes all 8.
    w=working(stats,rs,np.ones(8,bool));Q=[x['Q'] for x in w];q=[np.diag(x['Q'])*s['b'] for x,s in zip(w,stats)];v=[x['vp'] for x in w]
    f=custom(Q,q,v,rho=.8) if method=='Exact_rho08' else multisusie(Q,q,stats,v,L=5,rho=.75)
    f['pip']=np.where(keep,f['pip'],0.)
    assert np.isfinite(f['pip']).all() and np.isfinite(f['mean']).all()
    return f
def perturb(stats,rs,keep,random_seed):
    rng=np.random.default_rng(random_seed);w=working(stats,rs,keep);out=copy.deepcopy(stats)
    for k,c in enumerate(w):
        out[k]['b'][keep]+=np.linalg.cholesky(c['A'])@rng.normal(size=int(keep.sum()))/c['scale']
    return out
def jaccard(a,b):return float((a&b).sum()/max(1,(a|b).sum())) if (a|b).any() else 1.
def cell(phase,stage,region,pop):
    target=ROOT/'results'/phase/f'r{region}_{pop}';target.mkdir(parents=True,exist_ok=True)
    completion=target/f'{stage}_complete.json';stamp=codehash()
    if completion.exists():
        assert json.loads(completion.read_text())['code_sha256']==stamp,'Cache code mismatch'
        return f'{phase} {stage} r{region} {pop}: cached'
    cfg=CFG[phase];ls=lambdas(phase);arr,vs=load_reference(region)
    work=ROOT/'work'/phase/f'r{region}_{pop}'/stage;work.mkdir(parents=True,exist_ok=True)
    scenarios=['calibration'] if stage=='calibrate' else SCENARIOS
    if stage=='evaluate':
        thresholds=json.loads((target/'thresholds.json').read_text());assert thresholds['code_sha256']==stamp
        cuts=thresholds['cutoffs']
        if phase=='validation':assert (ROOT/'protocol/validation_threshold_freeze.json').exists()
    bank=[];cases=[]
    for scenario in scenarios:
        for rep in range(cfg[scenario]):
            stats,studyR,truth=generate(arr,pop,phase,region,scenario,rep)
            for ld in CFG['ld']:
                rs,fr=references(arr,pop,stats,studyR,ld);keep=eligible(stats,fr)
                id=f'{scenario}_{rep}_{ld}';sd=seed(phase,region,pop,scenario,rep)
                cases.append(case(stats,rs,keep,sd,id))
                item=dict(id=id,scenario=scenario,rep=rep,ld=ld,stats=stats,rs=rs,fr=fr,keep=keep,truth=truth,studyR=studyR)
                if phase=='validation' and stage=='evaluate' and rep<5:
                    cases.append(case(stats,rs,keep,sd+33000000,id+'_repeat'))
                    sp=perturb(stats,rs,keep,sd+44000000);item['perturbed']=sp
                    cases.append(case(sp,rs,keep,sd+55000000,id+'_perturb'))
                bank.append(item)
    start=time.perf_counter();lo=run_batch(cases,ls,work/'lasso');r_wall=time.perf_counter()-start
    cal=[];metrics=[];variants=[]
    for item in bank:
        id=item['id'];scenario=item['scenario'];rep=item['rep'];ld=item['ld'];s=item['stats'];rs=item['rs'];fr=item['fr'];keep=item['keep'];truth=item['truth'];causal=np.any(truth!=0,axis=0)
        score={};fits={}
        for l,lam in enumerate(ls):
            name=label(lam,phase);f=frequencies(lo[id],keep,l);score[name]=float(f.max())
            if stage=='evaluate':
                proposed=(f>cuts[ld][name])&keep
                fit=manc_run(work/'manc',s,vs,rs,fr,proposed,'joint');fit.update(frequency=f,proposed=proposed,converged=True)
                fit['seconds']+=lo[id]['seconds']/len(ls)
                fit['algorithmic_jaccard']='';fit['perturb_jaccard']='';fit['B200_jaccard']=''
                if 'perturbed' in item:
                    f2=frequencies(lo[id+'_repeat'],keep,l)
                    f3=frequencies(lo[id+'_perturb'],keep,l)
                    repeat=manc_run(work/'repeat',s,vs,rs,fr,(f2>cuts[ld][name])&keep,'joint')
                    pert=manc_run(work/'perturb',item['perturbed'],vs,rs,fr,(f3>cuts[ld][name])&keep,'joint')
                    fit['algorithmic_jaccard']=jaccard(fit['selected'],repeat['selected'])
                    fit['perturb_jaccard']=jaccard(fit['selected'],pert['selected'])
                    fit['B200_jaccard']=jaccard(proposed,((f+f2)/2>cuts[ld][name])&keep)
                fits[name]=fit
        score['Manc_stepwise']=entry_score(s,keep)
        if stage=='evaluate':
            f=manc_run(work/'baseline',s,vs,rs,fr,keep,'stepwise',cuts[ld]['Manc_stepwise'])
            f.update(proposed=f['selected'],converged=True,algorithmic_jaccard='',perturb_jaccard='',B200_jaccard='')
            if 'perturbed' in item:
                pert=manc_run(work/'baseline_perturb',item['perturbed'],vs,rs,fr,keep,'stepwise',cuts[ld]['Manc_stepwise'])
                f['perturb_jaccard']=jaccard(f['selected'],pert['selected'])
            fits['Manc_stepwise']=f
        for name in ['Exact_rho08','MultiSuSiE_L5']:
            st=time.perf_counter();f=bayesfit(s,rs,keep,name);elapsed=time.perf_counter()-st;score[name]=float(f['pip'].max())
            if not f['converged'] and stage=='calibrate':raise RuntimeError('Calibration Bayesian nonconvergence')
            if stage=='evaluate':
                sel=(f['pip']>cuts[ld][name])&keep
                fits[name]=dict(selected=sel,proposed=sel,mean=f['mean'],seconds=elapsed,status='ok' if f['converged'] else 'nonconverged',converged=f['converged'],dropped=0,algorithmic_jaccard='',perturb_jaccard='',B200_jaccard='')
        if stage=='calibrate':
            for name,value in score.items():cal.append(dict(ld=ld,method=name,replicate=rep,score=value))
            continue
        for name,f in fits.items():
            sel=f['selected'];n=int(sel.sum());nc=int(causal.sum());tp=int((sel&causal).sum());fp=int((sel&~causal).sum());err=f['mean']-truth
            metrics.append(dict(region=region,population=pop,scenario=scenario,replicate=rep,ld=ld,method=name,
              selected=n,proposed=int(f['proposed'].sum()),tp=tp,fp=fp,any_selection=int(n>0),fdp=fp/max(1,n),recall=tp/nc if nc else '',
              ld_recall_ceu=match_recovery(sel,causal,item['studyR'][0]) if nc else '',ld_recall_african=match_recovery(sel,causal,item['studyR'][1]) if nc else '',
              mse_ceu=float(np.mean(err[0]**2)),mse_african=float(np.mean(err[1]**2)),specific_estimate=float(f['mean'][1,6]) if scenario in ['heterogeneous','three'] else '',
              retained_causal=int((keep&causal).sum()),causal_count=nc,eligible_count=int(keep.sum()),joint_dropped=f['dropped'],
              seconds=f['seconds'],status=f['status'],converged=int(f['converged']),algorithmic_jaccard=f['algorithmic_jaccard'],perturb_jaccard=f['perturb_jaccard'],B200_jaccard=f['B200_jaccard']))
            for j in range(8):variants.append(dict(region=region,population=pop,scenario=scenario,replicate=rep,ld=ld,method=name,index=j,selected=int(sel[j]),truth_ceu=truth[0,j],truth_african=truth[1,j],estimate_ceu=f['mean'][0,j],estimate_african=f['mean'][1,j],frequency=f['frequency'][j] if 'frequency' in f else ''))
    if stage=='calibrate':
        table(target/'calibration_scores.tsv',cal);rank=math.ceil((cfg['calibration']+1)*.95)
        cut={ld:{name:float(np.sort([x['score'] for x in cal if x['ld']==ld and x['method']==name])[rank-1]) for name in score} for ld in CFG['ld']}
        save(target/'thresholds.json',dict(code_sha256=stamp,rank=rank,calibration_count=cfg['calibration'],lambdas=ls,cutoffs=cut))
    else:table(target/'metrics.tsv',metrics);table(target/'variants.tsv',variants)
    save(completion,dict(code_sha256=stamp,contexts=len(bank),R_batch_wall_seconds=r_wall,max_kkt=max(x['kkt'] for x in lo.values()),threshold_sha256=sha(target/'thresholds.json'),failed_fits=0))
    return f'{phase} {stage} r{region} {pop}: {len(bank)} contexts complete'
def freeze_tuning():
    import pandas as pd
    d=pd.concat([pd.read_csv(p,sep='\t',keep_default_na=False,na_values=['']) for p in (ROOT/'results/tuning').glob('*/metrics.tsv')]);d=d[d.ld=='matched'];choices=[]
    for lam in CFG['lambda_grid']:
        sub=d[d.method==label(lam,'tuning')];null=float(sub[sub.scenario=='null'].any_selection.mean());alt=sub[sub.scenario!='null']
        objective=float((alt.recall-alt.fdp).mean());assert np.isfinite(null) and np.isfinite(objective)
        choices.append(dict(lambda_value=lam,null_rate=null,recall_minus_fdp=objective))
    ok=[x for x in choices if x['null_rate']<=.10];chosen=max(ok,key=lambda x:(x['recall_minus_fdp'],x['lambda_value']))['lambda_value'] if ok else 2.5
    save(ROOT/'protocol/tuning_freeze.json',{'lambda':chosen,'candidates':choices,'fallback':not bool(ok),'code_sha256':codehash(),'evaluation_seed':CFG['validation']['seed'],'tuning_metrics_sha256':{str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'results/tuning').glob('*/metrics.tsv')}})
    print('Frozen lambda',chosen)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['tuning','validation'],required=True);ap.add_argument('--stage',choices=['calibrate','evaluate'],required=True);ap.add_argument('--workers',type=int,default=2);a=ap.parse_args()
    cfg=CFG[a.phase];regions=cfg.get('regions',CFG['regions']);pops=cfg.get('populations',CFG['populations'])
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futures=[ex.submit(cell,a.phase,a.stage,r,p) for r in regions for p in pops]
        for f in as_completed(futures):print(f.result(),flush=True)
    if a.phase=='tuning' and a.stage=='evaluate':freeze_tuning()
    if a.phase=='validation' and a.stage=='calibrate':
        save(ROOT/'protocol/validation_threshold_freeze.json',dict(code_sha256=codehash(),thresholds={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'results/validation').glob('*/thresholds.json')}))
if __name__=='__main__':main()
