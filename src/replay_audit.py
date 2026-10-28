"""Regenerate six prespecified final cases and replay both primary pipelines."""
from common import *
from lasso_bridge import case,run_batch,frequencies
from manc import run
from benchmark import seed
import pandas as pd
def main():
    arr,vs=load_reference(0);phase='validation';pop='YRI';region=0
    out=ROOT/'results/validation/r0_YRI';v=pd.read_csv(out/'variants.tsv',sep='\t',keep_default_na=False,na_values=[''])
    cuts=json.loads((out/'thresholds.json').read_text())['cutoffs'];lam=json.loads((ROOT/'protocol/tuning_freeze.json').read_text())['lambda']
    contexts=[];bank=[]
    for scenario in ['null','shared','heterogeneous']:
        stats,studyR,truth=generate(arr,pop,phase,region,scenario,0)
        for ld in CFG['ld']:
            rs,fr=references(arr,pop,stats,studyR,ld);keep=eligible(stats,fr);id=scenario+'_'+ld
            contexts.append(case(stats,rs,keep,seed(phase,region,pop,scenario,0),id));bank.append((id,scenario,ld,stats,rs,fr,keep))
    output=run_batch(contexts,[lam],ROOT/'work/replay/lasso');errors=[]
    for id,scenario,ld,stats,rs,fr,keep in bank:
        freq=frequencies(output[id],keep)
        for method in ['Bootstrap_LASSO_Manc','Manc_stepwise']:
            original=v[(v.scenario==scenario)&(v.replicate==0)&(v.ld==ld)&(v.method==method)].sort_values('index')
            if method.startswith('Bootstrap'):
                assert np.max(abs(original.frequency.astype(float).values-freq))<1e-12
                fit=run(ROOT/'work/replay/joint',stats,vs,rs,fr,keep&(freq>cuts[ld][method]),'joint')
            else:fit=run(ROOT/'work/replay/stepwise',stats,vs,rs,fr,keep,'stepwise',cuts[ld][method])
            assert np.array_equal(fit['selected'],original.selected.values.astype(bool))
            err=max(np.max(abs(fit['mean'][0]-original.estimate_ceu)),np.max(abs(fit['mean'][1]-original.estimate_african)))
            assert err<1e-9;errors.append(float(err))
    save(ROOT/'results/replay_audit.json',dict(passed=True,regenerated_contexts=6,primary_fits=12,frequency_exact_replay=True,max_coefficient_error=max(errors)))
    print('Six final contexts and both primary pipelines replayed successfully')
if __name__=='__main__':main()
