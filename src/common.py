"""Utilities; no personal R startup or environment files are read."""
from pathlib import Path
import os, sys, json, csv, hashlib
ROOT=Path(__file__).resolve().parents[1]
RUNTIME=os.environ.get('SCIENTIFIC_PYTHON','D:/Multiancestry_GWAS_Bayesian_Comparison/tmp/runtime')
sys.path.insert(0,RUNTIME)
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']: os.environ[k]='1'
import numpy as np
CFG=json.loads((ROOT/'protocol/design.json').read_text())
RSCRIPT=os.environ.get('RSCRIPT','C:/Program Files/R/R-4.5.1/bin/Rscript.exe')
def save(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x))+'\n',encoding='utf-8')
def table(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:return
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def codehash():
    files=[ROOT/'src'/p for p in ['common.py','bayes.py','manc.py','lasso_bridge.py','lasso.R','benchmark.py']]+[ROOT/'protocol/design.json']
    return hashlib.sha256(b''.join(p.read_bytes() for p in files)).hexdigest()
def load_reference(region):
    return dict(np.load(ROOT/f'data/local/region_{region}.npz')),json.loads((ROOT/f'data/local/region_{region}.json').read_text())
def generate(arr,pop,phase,region,scenario,rep,rows=False):
    sid=['calibration','null','shared','heterogeneous','one','three'].index(scenario)
    rng=np.random.default_rng(np.random.SeedSequence([CFG[phase]['seed'],region,CFG['populations'].index(pop),sid,rep]))
    truth=np.zeros((2,8));ix=CFG['causal_indices_zero_based'].get(scenario,[])
    if ix:truth[:,ix]=CFG['effects'][scenario]
    stats=[];rs=[];individuals=[]
    for k,(p,n) in enumerate(zip(['CEU',pop],CFG['cohort_sizes'])):
        h=arr[p+'_study'];i=rng.integers(len(h),size=(n,2));g=h[i[:,0]]+h[i[:,1]];x=g-g.mean(0)
        y=x@truth[k]+rng.normal(size=n);y-=y.mean();Q=x.T@x;d=np.diag(Q);q=x.T@y
        b=q/d;se=np.sqrt((y@y-q*q/d)/(n-2)/d)
        stats.append(dict(b=b,se=se,d=d,q=q,n=n,varY=float(y@y/(n-1)),freq=g.mean(0)/2))
        rs.append(Q/np.sqrt(d[:,None]*d[None,:]));individuals.append((g,y))
    return (stats,rs,truth,individuals) if rows else (stats,rs,truth)
def references(arr,pop,stats,studyR,ld):
    if ld=='in_sample':return studyR,np.array([s['freq'] for s in stats])
    return [arr[p+'_reference_R'] for p in ['CEU',pop]],np.array([arr[p+'_reference'].mean(0) for p in ['CEU',pop]])
def eligible(stats,fr):return np.all(np.abs(np.array([s['freq'] for s in stats])-fr)<CFG['frequency_limit'],axis=0)
def working(stats,rs,keep):
    """Mirror v1.1.0 GCTA scaling, including Vp from ALL input .ma lines."""
    out=[]
    for s,r in zip(stats,rs):
        v=2*s['freq']*(1-s['freq']);vp=float(np.median(v*s['n']*(s['se']**2+s['b']**2/(s['n']-1))))
        ne=(vp/v-s['b']**2)/s['se']**2+1;d=ne*v
        Q=r*np.minimum(ne[:,None],ne[None,:])*np.sqrt(v[:,None]*v[None,:]);np.fill_diagonal(Q,d)
        idx=np.flatnonzero(keep);H=Q[np.ix_(idx,idx)]/vp;scale=np.sqrt(np.diag(H));A=H/scale[:,None]/scale[None,:]
        z=scale*s['b'][idx]
        if len(idx) and np.linalg.eigvalsh(A).min()<=1e-8:raise ValueError('Non-PD working matrix; no silent repair')
        out.append(dict(A=A,z=z,scale=scale,Q=Q,vp=vp))
    return out
