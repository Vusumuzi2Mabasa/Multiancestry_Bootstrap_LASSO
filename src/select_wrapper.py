"""Standalone eight-marker research wrapper for harmonized Manc summary inputs.

This release requires two populations, eight common biallelic markers, complete
signed LD and a constant reported N within each population, matching validation.
No cutoff is guessed: users must supply an independently calibrated cutoff.
"""
from common import *
from lasso_bridge import case,run_batch,frequencies
from manc import run
import argparse
def readtab(path):
    lines=[x.split() for x in Path(path).read_text().splitlines() if x.strip()]
    return [dict(zip(lines[0],x)) for x in lines[1:]]
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sumstats',nargs=2,required=True);ap.add_argument('--ld',nargs=2,required=True)
    ap.add_argument('--lambda',dest='lam',type=float,required=True);ap.add_argument('--frequency-cutoff',type=float,required=True)
    ap.add_argument('--bootstraps',type=int,default=100);ap.add_argument('--seed',type=int,required=True);ap.add_argument('--out',required=True)
    a=ap.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    assert a.lam>0 and 0<=a.frequency_cutoff<=1 and a.bootstraps>0
    bim=[x.split() for x in Path(a.ld[0]+'.bim').read_text().splitlines() if x.strip()]
    assert len(bim)==8,'Validated research scope requires exactly eight BIM variants'
    vs=[dict(id=x[1],position=int(x[3]),alt=x[4],ref=x[5]) for x in bim];ids=[v['id'] for v in vs]
    assert len(set(ids))==8 and all(x[0]=='22' for x in bim),'This validated interface is chromosome 22 only'
    stats=[];rs=[];freq=[]
    for ma,pre in zip(a.sumstats,a.ld):
        rows=readtab(ma);lookup={x['SNP']:x for x in rows};assert set(lookup)==set(ids),'Provide the same eight markers to both cohorts'
        b=[];se=[];fr=[];ns=[]
        for v in vs:
            x=lookup[v['id']];sign=1 if (x['A1'],x['A2'])==(v['alt'],v['ref']) else -1
            assert (x['A1'],x['A2']) in [(v['alt'],v['ref']),(v['ref'],v['alt'])]
            f=float(x['freq']);fr.append(f if sign==1 else 1-f);b.append(sign*float(x['b']));se.append(float(x['se']));ns.append(float(x['N']))
        assert len(set(ns))==1 and ns[0]>2
        assert np.all(np.array(se)>0) and np.all((np.array(fr)>0)&(np.array(fr)<1))
        stats.append(dict(b=np.array(b),se=np.array(se),freq=np.array(fr),n=int(ns[0])))
        own=[x.split() for x in Path(pre+'.bim').read_text().splitlines() if x.strip()];orientation={x[1]:(x[4],x[5]) for x in own}
        sign=np.array([1 if orientation[v['id']]==(v['alt'],v['ref']) else -1 for v in vs])
        assert all(orientation[v['id']] in [(v['alt'],v['ref']),(v['ref'],v['alt'])] for v in vs)
        R=np.full((8,8),np.nan)
        for line in Path(pre+'.ld').read_text().splitlines():
            i,j,r=line.split();i=ids.index(i);j=ids.index(j);R[i,j]=R[j,i]=float(r)*sign[i]*sign[j]
        assert np.isfinite(R).all() and np.allclose(np.diag(R),1)
        rs.append(R);ff={x['SNP']:x for x in readtab(pre+'.frq')};af=[]
        for v in vs:
            x=ff[v['id']];f=float(x['MAF']);assert (x['A1'],x['A2']) in [(v['alt'],v['ref']),(v['ref'],v['alt'])]
            af.append(f if x['A1']==v['alt'] else 1-f)
        freq.append(af)
    fr=np.array(freq);keep=eligible(stats,fr)
    c=case(stats,rs,keep,a.seed,'input');lo=run_batch([c],[a.lam],out/'work',B=a.bootstraps)
    f=frequencies(lo['input'],keep);proposal=(f>a.frequency_cutoff)&keep
    result=run(out/'joint',stats,vs,rs,fr,proposal,'joint')
    table(out/'selection.tsv',[dict(variant=v['id'],eligible=int(keep[j]),bootstrap_frequency=f[j],proposed=int(proposal[j]),retained=int(result['selected'][j]),estimate_pop1=result['mean'][0,j],estimate_pop2=result['mean'][1,j]) for j,v in enumerate(vs)])
    save(out/'run.json',dict(method='Experimental bootstrap LASSO selection with official Manc joint refit',lambda_value=a.lam,frequency_cutoff=a.frequency_cutoff,B=a.bootstraps,seed=a.seed,kkt=lo['input']['kkt'],post_selection_intervals_valid=False,official_executable_modified=False))
    print(out/'selection.tsv')
if __name__=='__main__':main()
