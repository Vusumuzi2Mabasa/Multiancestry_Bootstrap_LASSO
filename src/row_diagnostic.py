from common import *
import subprocess
rows=[]
for region in [0,4]:
    arr,vs=load_reference(region)
    for scenario in ['null','shared','heterogeneous']:
        s,r,t,data=generate(arr,'YRI','tuning',region,scenario,980,True)
        rows.append(dict(region=region,scenario=scenario,rows=[dict(g=g,y=y) for g,y in data]))
save(ROOT/'work/row_input.json',rows)
r=subprocess.run([RSCRIPT,'--vanilla',str(ROOT/'src/row_diagnostic.R'),str(ROOT)],capture_output=True,text=True)
(ROOT/'work/row_diagnostic.log').write_text(r.stdout+'\n'+r.stderr)
if r.returncode:raise RuntimeError(r.stderr)
print('Genuine participant row-bootstrap diagnostic complete; private input remains ignored')
