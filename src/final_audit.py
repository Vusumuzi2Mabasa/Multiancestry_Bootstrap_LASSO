"""Reviewable publication manifest; no publishing or remote changes."""
from common import *
import subprocess,re,zipfile
def main():
    audit=json.loads((ROOT/'results/audit.json').read_text());assert audit['passed']
    assert json.loads((ROOT/'results/replay_audit.json').read_text())['passed']
    assert json.loads((ROOT/'results/cli_controls.json').read_text())['passed']
    readme=(ROOT/'README.md').read_text(encoding='utf-8')
    for text in ['Vusumuzi Mabasa','simulated','not an official','post-selection','FDR']:
        assert text.lower() in readme.lower()
    for target in re.findall(r'\]\(([^)]+)\)',readme):
        if not target.startswith(('http','#')):assert (ROOT/target).exists(),target
    doc=ROOT/'report/Bootstrap_LASSO_Manc_COJO_Research_Report.docx'
    with zipfile.ZipFile(doc) as z:
        core=z.read('docProps/core.xml').decode();xml=z.read('word/document.xml').decode()
        assert 'Vusumuzi Mabasa' in core and '<m:oMath' in xml
        assert 'NaN' not in xml and 'nan (' not in xml
    assert (ROOT/'work/render_github/page-1.png').exists(),'Report must be rendered before final audit'
    visual=json.loads((ROOT/'work/visual_review_github.json').read_text())
    assert visual['passed'] and visual['pages_inspected']==len(list((ROOT/'work/render_github').glob('page-*.png')))
    tracked=subprocess.check_output(['git','-C',str(ROOT),'ls-files'],text=True).splitlines()
    assert not any(p.startswith(('data/local/','work/')) or p.endswith(('.npz','.bed','.fam','.vcf','.vcf.gz','.exe')) for p in tracked)
    assert not any(Path(p).name in ['.Renviron','.Rprofile','.Rhistory','.RData'] for p in tracked)
    for name in tracked:
        if Path(name).suffix.lower() in ['.py','.r','.md','.txt','.tsv','.json']:
            content=(ROOT/name).read_text(encoding='utf-8')
            assert not re.search(r'\b(?:HG\d{5}|NA\d{5})\b',content),'Possible donor identifier in '+name
    remote_urls=subprocess.check_output(['git','-C',str(ROOT),'remote','-v'],text=True).splitlines()
    expected_remote='https://github.com/Vusumuzi2Mabasa/Multiancestry_Bootstrap_LASSO.git'
    assert all(line.split()[1]==expected_remote for line in remote_urls),'Unexpected publication destination'
    predecessor=Path('D:/Multiancestry_GWAS_Bayesian_Comparison')
    assert subprocess.check_output(['git','-C',str(predecessor),'rev-parse','HEAD'],text=True).strip()==CFG['source_study_commit']
    assert not subprocess.check_output(['git','-C',str(predecessor),'status','--porcelain'],text=True).strip()
    manifest=[]
    for name in sorted(tracked):
        if name not in ['checksums_public.sha256','results/final_review.json']:
            manifest.append(f'{sha(ROOT/name)}  {name}')
    (ROOT/'checksums_public.sha256').write_text('\n'.join(manifest)+'\n',encoding='utf-8',newline='\n')
    save(ROOT/'results/final_review.json',dict(passed=True,raw_genotypes_tracked=False,donor_identifier_files_tracked=False,remote_configured=bool(remote_urls),
      report_author='Vusumuzi Mabasa',native_equations=True,all_results_audited=True,manifest_files=len(manifest),
      source_study_preserved_commit=CFG['source_study_commit'],visual_review='All rendered report pages inspected; see work/visual_review_github.json'))
    print('Final publication audit passed')
if __name__=='__main__':main()
