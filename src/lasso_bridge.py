from common import *
import subprocess
def run_batch(cases,lambdas,work,B=None):
    work=Path(work);work.mkdir(parents=True,exist_ok=True)
    src=work/'input.json';dst=work/'output.json'
    save(src,dict(B=B or CFG['bootstrap_B'],lambdas=lambdas,cases=cases))
    r=subprocess.run([RSCRIPT,'--vanilla',str(ROOT/'src/lasso.R'),str(src),str(dst)],capture_output=True,text=True)
    (work/'R.log').write_text(r.stdout+'\n'+r.stderr)
    if r.returncode:raise RuntimeError(r.stderr[-2500:])
    return {c['id']:c for c in json.loads(dst.read_text())}
def case(stats,rs,keep,seed,id):
    w=working(stats,rs,keep)
    return dict(id=id,seed=int(seed),context=[dict(A=c['A'],z=c['z']) for c in w])
def frequencies(output,keep,column=0):
    f=np.zeros(8);p=int(keep.sum())
    if p:f[keep]=np.array(output['frequencies']).reshape(p,-1)[:,column]
    return f
