"""Reuse local public-reference derivatives, verify hashes, never copy donor IDs."""
from common import *
import argparse,shutil,urllib.request
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-study',required=True);a=ap.parse_args();source=Path(a.source_study)
    provenance=json.loads((ROOT/'results/provenance.json').read_text())
    (ROOT/'data/local').mkdir(parents=True,exist_ok=True)
    for record in provenance['reference_files']:
        dest=ROOT/record['file'];src=source/'data/revision/reference'/dest.name
        assert sha(src)==record['sha256'],f'Input identity mismatch: {dest.name}'
        if dest.exists():assert sha(dest)==record['sha256']
        else:shutil.copyfile(src,dest)
    exe=ROOT/'tools/manc_cojo_win.exe';local=source/'tools/manc_cojo_win.exe'
    if not exe.exists():
        if local.exists():shutil.copyfile(local,exe)
        else:urllib.request.urlretrieve(provenance['manc_url'],exe)
    assert sha(exe)==provenance['manc_sha256']
    print('Verified six local public-reference panels and official Manc executable')
if __name__=='__main__':main()
