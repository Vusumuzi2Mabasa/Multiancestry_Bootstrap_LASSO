"""Unmodified official executable, common inputs, selection or joint-only mode."""
from common import *
import subprocess, time, math
from scipy.special import ndtr
def entry_score(stats,keep):
    if not keep.any():return 0.
    w=np.array([1/s['se']**2 for s in stats]);b=np.array([s['b'] for s in stats])
    z=(w*b).sum(0)/np.sqrt(w.sum(0))
    return -math.log10(max(float((2*ndtr(-abs(z[keep]))).min()),1e-300))
def write_inputs(work,stats,variants,rs,fr):
    work=Path(work);work.mkdir(parents=True,exist_ok=True);refs=[];mas=[]
    for k in range(2):
        pre=work/f'ld{k}';refs.append(str(pre));ma=work/f's{k}.ma';mas.append(str(ma))
        pre.with_suffix('.bim').write_text(''.join(f"22\t{v['id']}\t0\t{v['position']}\t{v['alt']}\t{v['ref']}\n" for v in variants))
        pre.with_suffix('.ld').write_text(''.join(f"{a['id']} {b['id']} {rs[k][i,j]:.17g}\n" for i,a in enumerate(variants) for j,b in enumerate(variants) if j>=i))
        table(pre.with_suffix('.frq'),[dict(CHR=22,SNP=v['id'],A1=v['alt'],A2=v['ref'],MAF=fr[k,j],NCHROBS=80) for j,v in enumerate(variants)])
        s=stats[k]
        # ALL marginal rows retained to keep Vp unchanged between select and joint.
        table(ma,[dict(SNP=v['id'],A1=v['alt'],A2=v['ref'],freq=s['freq'][j],b=s['b'][j],se=s['se'][j],p=max(float(2*ndtr(-abs(s['b'][j]/s['se'][j]))),1e-300),N=s['n']) for j,v in enumerate(variants)])
    return refs,mas
def run(work,stats,variants,rs,fr,keep,mode,cutoff=None,screen=True):
    start=time.perf_counter();sel=np.zeros(8,bool);eff=np.zeros((2,8));se=np.zeros((2,8))
    if not keep.any():return dict(selected=sel,mean=eff,se=se,status='empty',seconds=time.perf_counter()-start,dropped=0)
    if mode=='stepwise' and screen and entry_score(stats,keep)<=cutoff:
        return dict(selected=sel,mean=eff,se=se,status='analytic_no_entry',seconds=time.perf_counter()-start,dropped=0)
    refs,mas=write_inputs(work,stats,variants,rs,fr);work=Path(work)
    extract=work/'extract.txt';extract.write_text(''.join(v['id']+'\n' for j,v in enumerate(variants) if keep[j]))
    output=work/'fit.jma.cojo'
    if output.exists():output.unlink()
    cmd=[str(ROOT/'tools/manc_cojo_win.exe'),'--ld',*refs,'--cojo-file',*mas,'--extract',str(extract),'--diff-freq','1','1','--thread-num','1','--effect-mode','GCTA','--out',str(work/'fit')]
    cmd+=['--cojo-slct','--cojo-p',str(10**(-cutoff))] if mode=='stepwise' else ['--cojo-joint']
    result=subprocess.run(cmd,capture_output=True,text=True,timeout=45)
    log=result.stdout+'\n'+result.stderr;(work/'last_stdout.txt').write_text(log)
    (work/'command.json').write_text(json.dumps(cmd))
    if result.returncode:raise RuntimeError(log[-2000:])
    if output.exists():
        lines=[x.split() for x in output.read_text().splitlines() if x.strip()];lookup={v['id']:j for j,v in enumerate(variants)}
        for row in [dict(zip(lines[0],x)) for x in lines[1:]]:
            j=lookup[row['SNP']];assert row['A1']==variants[j]['alt'];sel[j]=True
            eff[:,j]=[float(row['bJ.1']),float(row['bJ.2'])];se[:,j]=[float(row['bJ_se.1']),float(row['bJ_se.2'])]
    elif not any(x in log.lower() for x in ['no snp','no significant','0 snp','no variant','not significant','no joint snp']):raise RuntimeError(log[-1500:])
    assert np.isfinite(eff).all() and not np.any(sel&~keep)
    return dict(selected=sel,mean=eff,se=se,status='official_'+mode,seconds=time.perf_counter()-start,dropped=int(keep.sum()-sel.sum()) if mode=='joint' else 0)
