from common import *
from manc import run,entry_score
from lasso_bridge import run_batch,case,frequencies
import subprocess,time
def main():
    r=subprocess.run([RSCRIPT,'--vanilla',str(ROOT/'src/test_lasso.R'),str(ROOT)],capture_output=True,text=True)
    (ROOT/'work/test_R.log').write_text(r.stdout+'\n'+r.stderr)
    if r.returncode:raise RuntimeError(r.stderr)
    arr,vs=load_reference(0);s,rs,truth=generate(arr,'YRI','tuning',0,'shared',900)
    fr=np.array([x['freq'] for x in s]);keep=np.ones(8,bool);w=working(s,rs,keep)
    errs=[];seerrs=[];replays=[]
    for ix in [[1],[1,6],[0,2,5]]:
        use=np.isin(np.arange(8),ix)
        f=run(ROOT/'work/tests/joint',s,vs,rs,fr,use,'joint')
        actual=np.flatnonzero(f['selected'])
        for k in range(2):
            Q=w[k]['Q'][np.ix_(actual,actual)];q=np.diag(w[k]['Q'])[actual]*s[k]['b'][actual]
            beta=np.linalg.solve(Q,q);sd=np.sqrt(np.diag(np.linalg.inv(Q))*w[k]['vp'])
            if len(actual)==1:sd=s[k]['se'][actual]
            errs.extend(abs(f['mean'][k,actual]-beta));seerrs.extend(abs(f['se'][k,actual]-sd))
    for scenario in ['null','shared','heterogeneous']:
        s,rs,truth=generate(arr,'YRI','tuning',0,scenario,901);fr=np.array([x['freq'] for x in s])
        f=run(ROOT/'work/tests/step',s,vs,rs,fr,keep,'stepwise',1.3,False)
        j=run(ROOT/'work/tests/replay',s,vs,rs,fr,f['selected'],'joint')
        replays.append(float(np.max(abs(f['mean']-j['mean']))))
        screened=run(ROOT/'work/tests/screen',s,vs,rs,fr,keep,'stepwise',10)
        unscreened=run(ROOT/'work/tests/unscreen',s,vs,rs,fr,keep,'stepwise',10,False)
        assert np.array_equal(screened['selected'],unscreened['selected'])
    assert max(errs)<1e-5 and max(seerrs)<1e-5 and max(replays)<1e-5
    c=case(s,rs,keep,123,'a');o=run_batch([c],[1.5],ROOT/'work/tests/batch')
    assert np.all((frequencies(o['a'],keep)>=0)&(frequencies(o['a'],keep)<=1))
    save(ROOT/'results/pipeline_controls.json',dict(passed=True,joint_max_error=max(errs),joint_se_max_error=max(seerrs),stepwise_joint_replay_max_error=max(replays),first_entry_screen_passed=True))
    print('R solver and official Manc joint/selection controls passed')
if __name__=='__main__':main()
