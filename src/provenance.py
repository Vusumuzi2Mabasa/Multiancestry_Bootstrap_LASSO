from common import *
import shutil,subprocess,platform
def main():
    source=Path('D:/Multiancestry_GWAS_Bayesian_Comparison')
    files=[]
    for p in sorted((ROOT/'data/local').glob('region_*')):
        original=source/'data/revision/reference'/p.name
        assert sha(p)==sha(original)
        files.append(dict(file='data/local/'+p.name,sha256=sha(p)))
    save(ROOT/'results/provenance.json',dict(source_study_commit=CFG['source_study_commit'],upstream_manc_commit=CFG['manc_source_commit'],
      reference_files=files,manc_sha256=sha(ROOT/'tools/manc_cojo_win.exe'),
      public_vcf_sha256='a90c16c4ff2b3196476d506ae13cb3047fae8670163c7c932c4b0239aef3daf5',
      source_study='https://github.com/Vusumuzi2Mabasa/Multiancestry_GWAS_Bayesian_Comparison',
      manc_url='https://github.com/light156/multi-ancestry-COJO/releases/download/v1.1.0/manc_cojo_win.exe',
      population_donors={'CEU':{'study':59,'reference':40},'YRI':{'study':68,'reference':40},'LWK':{'study':59,'reference':40}},
      raw_data_published=False,phenotypes='SIMULATED',python=platform.python_version(),R_startup='--vanilla',
      numpy=np.__version__,multisusie_commit='a89e091d45297d983118438f116c3eb5bc1219ed'))
if __name__=='__main__':main()
