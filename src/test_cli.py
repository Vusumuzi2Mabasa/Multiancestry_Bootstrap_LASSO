"""Allele-flip invariance of the user-facing Manc-format wrapper."""
from common import *
import subprocess,csv
def main():
    source=ROOT/'work/tests/joint';flip=ROOT/'work/cli_flip';flip.mkdir(parents=True,exist_ok=True)
    for k in range(2):
        rows=list(csv.DictReader((source/f's{k}.ma').open(),delimiter='\t'));marker=rows[0]['SNP']
        for row in rows:
            if row['SNP']==marker:
                row['A1'],row['A2']=row['A2'],row['A1'];row['b']=-float(row['b']);row['freq']=1-float(row['freq'])
        table(flip/f's{k}.ma',rows)
        bim=[l.split() for l in (source/f'ld{k}.bim').read_text().splitlines()]
        for row in bim:
            if row[1]==marker:row[4],row[5]=row[5],row[4]
        (flip/f'ld{k}.bim').write_text('\n'.join('\t'.join(row) for row in bim)+'\n')
        ld=[l.split() for l in (source/f'ld{k}.ld').read_text().splitlines()]
        for row in ld:
            if (row[0]==marker)!=(row[1]==marker):row[2]=str(-float(row[2]))
        (flip/f'ld{k}.ld').write_text('\n'.join(' '.join(row) for row in ld)+'\n')
        fr=list(csv.DictReader((source/f'ld{k}.frq').open(),delimiter='\t'))
        for row in fr:
            if row['SNP']==marker:row['A1'],row['A2']=row['A2'],row['A1'];row['MAF']=1-float(row['MAF'])
        table(flip/f'ld{k}.frq',fr)
    # Parametric random draws are distributionally invariant but not identical under
    # a Cholesky allele flip with the same random normals. Compare near-unpenalized
    # full support and joint estimates to isolate deterministic file harmonization.
    outputs=[]
    for path,label in [(source,'original'),(flip,'flipped')]:
        out=ROOT/'work/cli_control'/label
        cmd=[sys.executable,str(ROOT/'src/select_wrapper.py'),'--sumstats',str(path/'s0.ma'),str(path/'s1.ma'),'--ld',str(path/'ld0'),str(path/'ld1'),'--lambda','.000001','--frequency-cutoff','0','--seed','731593','--out',str(out)]
        r=subprocess.run(cmd,capture_output=True,text=True);assert r.returncode==0,r.stderr
        outputs.append(list(csv.DictReader((out/'selection.tsv').open(),delimiter='\t')))
    errors=[]
    for j,(a,b) in enumerate(zip(*outputs)):
        assert a['retained']==b['retained'] and a['eligible']==b['eligible']
        for col in ['estimate_pop1','estimate_pop2']:
            errors.append(abs(float(a[col])*(-1 if j==0 else 1)-float(b[col])))
    assert max(errors)<1e-7
    save(ROOT/'results/cli_controls.json',dict(passed=True,control='Full-support allele flip with deterministic joint-estimator comparison',max_effect_error=max(errors),bootstrap_draw_by_draw_invariance_claimed=False))
    print('Manc-format wrapper allele orientation control passed')
if __name__=='__main__':main()
