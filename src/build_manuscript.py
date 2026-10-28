"""Generate the research README from audited results, with interpretation first."""
from common import *
import pandas as pd
def f(x,n=3):return f'{float(x):.{n}f}'
def main():
    audit=json.loads((ROOT/'results/audit.json').read_text());assert audit['passed']
    s=pd.read_csv(ROOT/'results/validation_summary.tsv',sep='\t',keep_default_na=False,na_values=['']);p=pd.read_csv(ROOT/'results/paired_differences.tsv',sep='\t',keep_default_na=False,na_values=[''])
    t=pd.read_csv(ROOT/'results/threshold_summary.tsv',sep='\t');tune=json.loads((ROOT/'protocol/tuning_freeze.json').read_text())
    def r(method,sc,ld='matched'):return s[(s.method==method)&(s.scenario==sc)&(s.ld==ld)].iloc[0]
    def difference(sc,metric,ld='matched'):
        a=p[(p.comparator=='Manc_stepwise')&(p.scenario==sc)&(p.metric==metric)&(p.ld==ld)].iloc[0]
        return f"{f(a.difference)} ({f(a.lower)} to {f(a.upper)})"
    boot='Bootstrap_LASSO_Manc';manc='Manc_stepwise';names={boot:'Bootstrap LASSO plus Manc joint',manc:'Original Manc stepwise','Exact_rho08':'Exact Bayes rho 0.8','MultiSuSiE_L5':'MultiSuSiE L5'}
    nulltable='| Method | Null rate | MC SE | Wilson interval |\n| --- | --- | --- | --- |\n'
    for m,n in names.items():
        x=r(m,'null');nulltable+=f'| {n} | {f(x.any_selection)} | {f(x.any_selection_mcse)} | {f(x.null_wilson_lower)} to {f(x.null_wilson_upper)} |\n'
    perf='| Scenario and method | Exact recall | LD recall AFR | Mean FDP | AFR RMSE |\n| --- | --- | --- | --- | --- |\n'
    for sc in ['one','shared','heterogeneous','three']:
        for m,n in [(boot,'Bootstrap'),(manc,'Stepwise')]:
            x=r(m,sc);perf+=f'| {sc.capitalize()} {n} | {f(x.recall)} | {f(x.ld_recall_african)} | {f(x.fdp)} | {f(x.rmse_african)} |\n'
    paired='| Scenario | Recall change and 95% MC interval | FDP change and 95% MC interval |\n| --- | --- | --- |\n'
    for sc in ['one','shared','heterogeneous','three']:paired+=f'| {sc.capitalize()} | {difference(sc,"recall")} | {difference(sc,"fdp")} |\n'
    tuning='| Lambda | Tuning null rate | Recall minus FDP |\n| --- | --- | --- |\n'
    for x in tune['candidates']:tuning+=f'| {x["lambda_value"]} | {f(x["null_rate"])} | {f(x["recall_minus_fdp"])} |\n'
    stable='| Scenario | Bootstrap repeat Jaccard | Bootstrap perturb Jaccard | Stepwise perturb Jaccard |\n| --- | --- | --- | --- |\n'
    for sc in ['null','one','shared','heterogeneous','three']:
        a=r(boot,sc);b=r(manc,sc);stable+=f'| {sc.capitalize()} | {f(a.algorithmic_jaccard)} | {f(a.perturb_jaccard)} | {f(b.perturb_jaccard)} |\n'
    cb=t[(t.method==boot)&(t.ld=='matched')].cutoff
    row=pd.read_csv(ROOT/'results/row_bootstrap_diagnostic.tsv',sep='\t');rd=float(abs(row.row_frequency-row.parametric_frequency).mean());rm=float(abs(row.row_frequency-row.parametric_frequency).max())
    tests=json.loads((ROOT/'results/lasso_controls.json').read_text());joint=json.loads((ROOT/'results/pipeline_controls.json').read_text())
    methods=(ROOT/'protocol/methods.md').read_text().split('## Retained interfaces and changes',1)[1]
    methods='### Retained interfaces and changes'+methods
    methods=methods.replace('\n## ','\n### ')
    text=f'''# Replacing Manc COJO stepwise selection with bootstrap LASSO

**Vusumuzi Mabasa**

This study tests an experimental replacement for the selection stage of Manc-COJO. It uses public CEU, YRI and LWK reference haplotypes to generate study genotypes; **every phenotype is simulated**. No AWI-Gen data, observed clinical outcomes or biological discoveries are claimed. The new pipeline is a wrapper using R bootstrap LASSO selection followed by the unchanged official Manc-COJO joint estimator. It is not an official release, forked executable or endorsement by the Manc-COJO authors.

[Word manuscript](report/Bootstrap_LASSO_Manc_COJO_Research_Report.docx) · [Protocol](protocol/design.json) · [Mathematical specification](protocol/methods.md) · [Audited results](results/validation_summary.tsv) · [Paired uncertainty](results/paired_differences.tsv) · [Workflow](workflow.txt)

## Abstract

**Background.** Replacing stepwise signal selection with bootstrap LASSO changes the statistical selection rule, even when summary inputs and subsequent joint estimation are held fixed. The relevant comparison is therefore between complete pipelines with a shared estimation interface.

**Methods.** The primary comparison used original Manc-COJO v1.1.0 selection versus an R ordinary LASSO selector fitted to a parametric bootstrap of ancestry-specific summary likelihoods. Both used the same harmonized marginal statistics, signed LD and candidate variants, and the official joint estimator. Independent tuning selected lambda {tune['lambda']}; separate null simulations calibrated selection cutoffs before fresh validation. Validation comprised 2,388 calibration and 2,160 evaluation datasets across 12 fixed strata, analysed under matched external and in-sample LD. Exact Bayesian and MultiSuSiE models were secondary comparators.

**Results.** With matched external LD, bootstrap selection had a held-out global-null rate of {f(r(boot,'null').any_selection)}, compared with {f(r(manc,'null').any_selection)} for original stepwise selection. Shared-effect exact recall and mean FDP were {f(r(boot,'shared').recall)} and {f(r(boot,'shared').fdp)} for the replacement, versus {f(r(manc,'shared').recall)} and {f(r(manc,'shared').fdp)} for the original. In the heterogeneous setting, the paired recall change was {difference('heterogeneous','recall')}. All {audit['paired_method_rows']:,} evaluation method rows passed completeness and numerical reconstruction checks.

**Conclusion.** The replacement improved mean recovery in heterogeneous and three-causal settings, but reduced recovery in one-causal and shared-effect settings and required more computation. It did not consistently improve perturbation stability. Bootstrap frequencies are not posterior inclusion probabilities, valid post-selection confidence intervals or a guarantee of nonnull FDR control.

## Introduction

Conditional and joint GWAS analysis aims to separate association signals in the presence of linkage disequilibrium. Manc-COJO extends this analysis across populations with different LD structures [1]. Its selection stage uses conditional association evidence, while the final joint stage estimates effects for a retained set. A change in the selector should be assessed while preserving the surrounding data handling and estimation steps wherever feasible.

The research question is whether replacing Manc-COJO's forward and backward selection with bootstrap LASSO gives useful differences in error, recovery and stability. Ordinary LASSO produces sparse coefficient vectors; repeated bootstrap fits provide an empirical selection frequency [2,3]. That frequency describes sensitivity to the stated resampling model. It is not interchangeable with a Bayesian probability or an error-controlled stability-selection procedure.

This project is independent of the [preceding public-reference study](https://github.com/Vusumuzi2Mabasa/Multiancestry_GWAS_Bayesian_Comparison/tree/03a5a07826a82380e6a1a9191d4f0a7f61f0f6af). It reuses its disjoint reference panels and simulated-effect architectures but uses new seeds and a new frozen protocol. The original project remains intact. All conclusions concern quantitative simulated traits in six sparse chromosome 22 regions. YRI and LWK identify specific source panels, not African genetic diversity as a whole.

## Methods

{methods}

## Results

### Numerical agreement and independent diagnostics

The R coordinate-descent solver agreed with glmnet to a maximum absolute coefficient difference of {tests['glmnet_max_difference']:.2g}. Orthogonal soft-threshold, allele-flip, marker-permutation and KKT controls passed. A 50,000-draw covariance control differed from its target by at most {f(tests['covariance_max_difference'])}; the independent fixed-design residual simulation also passed its prespecified covariance tolerance. Production KKT residuals were at most {audit['max_kkt']:.2g}.

The independently reconstructed GCTA joint coefficients differed from the official executable by at most {joint['joint_max_error']:.2g}; joint SE agreement was {joint['joint_se_max_error']:.2g}. Replaying stepwise selected sets in joint-only mode gave maximum coefficient difference {joint['stepwise_joint_replay_max_error']:.2g}. These controls validate retained calculations; they do not validate post-selection inference or external-LD assumptions.

Six separate participant-level diagnostic datasets compared 500 genuine within-population row bootstraps with 500 summary parametric draws at lambda 1.5 and known noise variance one. Across 48 markers, mean absolute frequency difference was {f(rd)} and the maximum was {f(rm)}. This is a descriptive diagnostic on a different, oracle-variance specification, not a substitution for the primary summary-input experiment or proof that the two bootstraps are equivalent.

### Independent tuning and calibration

The prespecified tuning rule selected lambda {tune['lambda']}. Development results below use only two tuning regions paired with YRI. They were not pooled with validation, and final evaluation outcomes did not affect the choice. The lowest penalty can saturate bootstrap frequencies and yield a cutoff of one, causing an empty selector; this negative behavior was retained rather than bypassed.

{tuning}

Table 1. Tuning results under matched external LD. Recall minus FDP is the prespecified development objective, not a clinical utility score or an FDR guarantee.

Fresh matched-LD bootstrap cutoffs ranged from {f(cb.min())} to {f(cb.max())} across the 12 strata. Every cutoff was saved and hashed before validation evaluation. Held-out null results use 720 datasets per method; rates estimate locus-level any-selection error conditional on the fixed panels and calibration thresholds. Pooled Wilson intervals are descriptive approximations; the primary MC SE uses the fixed-stratum formula and allows different rates across strata.

{nulltable}

Table 2. Fresh matched-reference global-null error. No genome-wide multiplicity or nonnull FDR interpretation is intended.

### Paired recovery and false discoveries

The shared and heterogeneous scenarios distinguish recovery of broadly shared effects from recovery of a population-specific effect in the smaller cohort. For shared effects, bootstrap recall changed by {difference('shared','recall')} relative to stepwise selection, and FDP changed by {difference('shared','fdp')}. In the heterogeneous setting, the corresponding FDP change was {difference('heterogeneous','fdp')}. These paired intervals use differences within the same simulated dataset and then aggregate equally across fixed strata.

{perf}

Table 3. Primary matched-LD results, 360 datasets per alternative and method. AFR is the African-source cohort. RMSE covers all eight coefficients, assigning zero to unselected markers. LD-aware recovery uses one-to-one matching and the original causal denominator.

{paired}

Table 4. Replacement minus original stepwise selection. Positive recall changes indicate more causal retention; positive FDP changes indicate a larger false discovery proportion. Intervals quantify Monte Carlo uncertainty conditional on frozen settings.

The in-sample LD diagnostic changes the information supplied to both pipelines. Shared-effect FDP was {f(r(boot,'shared','in_sample').fdp)} for the replacement and {f(r(manc,'shared','in_sample').fdp)} for the original, compared with {f(r(boot,'shared').fdp)} and {f(r(manc,'shared').fdp)} under matched external LD. Thus error patterns must be interpreted jointly with LD source; this comparison does not show that bootstrap perturbation accounts for reference-panel uncertainty.

![Primary paired comparison with simulated phenotypes](figures/primary_comparison.png)

Figure 1. Primary null error, exact recovery and mean FDP. Solid curves and circles use matched external LD; dashed curves and lighter squares use in-sample LD. Error bars are 1.96 Monte Carlo SEs. Global-null calibration does not constrain alternative FDP.

### Effect estimation and secondary Bayesian comparisons

Both primary pipelines use the official joint estimator on their own selected sets. The mean estimated African-source effect at the population-specific marker, whose generating effect is 0.28, was {f(r(boot,'heterogeneous').specific_estimate)} for bootstrap selection and {f(r(manc,'heterogeneous').specific_estimate)} for stepwise selection. These averages include zero when the marker was not retained; they combine selection loss and estimation error. The paired African-source all-marker MSE difference in this scenario was {difference('heterogeneous','mse_african')}. There are no post-selection confidence-interval claims.

Exact Bayes and MultiSuSiE provide secondary comparisons with different priors and estimators. In the shared matched-LD setting, their recall was {f(r('Exact_rho08','shared').recall)} and {f(r('MultiSuSiE_L5','shared').recall)}, with mean FDP {f(r('Exact_rho08','shared').fdp)} and {f(r('MultiSuSiE_L5','shared').fdp)}, respectively. Their posterior means use shrinkage across all markers, whereas the two primary pipelines refit selected sets. These comparisons contextualize the tradeoffs but do not isolate the selector as directly as the original-versus-replacement comparison.

![Secondary methods with simulated phenotypes](figures/secondary_comparison.png)

Figure 2. Secondary matched-reference comparisons of recovery, FDP and coefficient RMSE. Bayesian PIPs and bootstrap frequencies have distinct meanings and independently calibrated cutoffs. Full numerical estimates and paired uncertainty are supplied in the result tables.

### Stability runtime and fit accounting

Bootstrap repeat stability measures Monte Carlo sensitivity of the new selector on unchanged data. Common perturbation stability applies the same marginal-statistic perturbation to both pipelines. Each stability entry averages 60 paired cases per scenario and LD condition, selected by replicate index before observing results. A null empty-empty comparison scores one, so null Jaccard values should not be interpreted as evidence of informative recovery.

{stable}

Table 5. Matched-LD stability. Algorithmic repeat Jaccard and common perturbation Jaccard answer different questions. The B200 diagnostic is available in the complete summary and uses the original B100 cutoff without recalibration.

Bootstrap aggregation did not automatically improve stability. In heterogeneous matched-LD cases, common-perturbation Jaccard was {f(r(boot,'heterogeneous').perturb_jaccard)} for the replacement versus {f(r(manc,'heterogeneous').perturb_jaccard)} for stepwise selection. Greater causal recovery therefore coexisted with lower stability under this particular perturbation. The perturbation experiment is not a held-out global-null test and has no claimed error calibration.

For shared matched-LD cases, measured average runtime was {f(r(boot,'shared').seconds)} seconds for the replacement and {f(r(manc,'shared').seconds)} seconds for original stepwise selection. Replacement timings include R fitting time and official joint-process wall time; batched R startup and serialization are recorded separately in shard completion files. Repeated stability fits are excluded from primary timings. Null original fits may use a numerically checked no-entry screen. These implementation timings include file and process overhead, depend on the Toshiba drive and concurrent workloads, and are not a universal algorithmic speed ranking.

![Stability and runtime with simulated phenotypes](figures/stability_runtime.png)

Figure 3. Common-perturbation stability and measured runtime under matched external LD. Both pipelines receive the same perturbed summary input in each stability comparison.

The audit retained all {audit['paired_method_rows']:,} evaluation rows and {audit['variant_rows']:,} per-variant records. There were {audit['failed_fits']} failed fits and {audit['nonconverged']} finite nonconverged Bayesian fits. Official joint safeguards removed {audit['total_joint_drops']} proposed markers across the recorded primary replacement contexts. Causal-filter accounting and proposed-versus-final counts remain in each metric row; recall never drops filtered causal variants from its denominator. Reconstructed MSEs agreed with saved metrics to {audit['max_mse_reconstruction_error']:.2g}.

## Discussion

The clearest favorable recovery changes occurred in heterogeneous and three-causal settings, with paired recall increases of {difference('heterogeneous','recall')} and {difference('three','recall')}. The one-causal setting gave an unfavorable combination: recall changed by {difference('one','recall')} and FDP by {difference('one','fdp')}. Shared-effect recall also decreased. The heterogeneous FDP interval spans zero, so its positive mean difference alone does not establish a reproducible increase at this simulation precision. These findings argue against replacing the original selector unconditionally.

Replacing stepwise selection is a change in the objective and in the way evidence is combined, not merely a different search for an identical optimum. The specified LASSO uses separate population coefficients and selects the union of their supports. It can retain a population-specific signal without strong fixed-effect meta-analysis evidence, but it also creates additional routes for noise to enter. Independent global-null calibration makes the comparison more interpretable without forcing the same alternative FDP.

The parametric bootstrap makes the summary-input replacement possible, but it conditions on the supplied LD and estimated variance. It cannot repair biased external LD by repeatedly drawing from that same working model. The saturated pilot also propagates observed noise. The in-sample control and participant-level diagnostic help locate these assumptions; neither turns the primary frequencies into calibrated posterior probabilities.

The preserved joint estimator is a strength of the comparison. Any change in the primary final estimates arises from the selected set and the official safeguards, not an undisclosed switch to a new coefficient estimator. Its unadjusted SEs remain unsuitable as a general post-selection inferential guarantee. A scientifically useful result may therefore be a recovery improvement accompanied by extra FDP, or a failure of bootstrap selection to improve a target metric. The tables retain both directions of change.

## Limitations

Only eight variants per region and six chromosome 22 windows are analysed. The fixed reference pools are small, study participants are resampled haplotype pairs, and there is no admixture, relatedness, binary trait, genotype missingness or overlapping-cohort model. Sparse positive-definite panels permit an unregularized saturated pilot; dense or singular regions require a newly specified method and validation. The standalone interface is intentionally restricted to the tested two-cohort, eight-marker, constant-N setting.

Lambda tuning covered two regions and one African-source panel. Fresh seeds prevent outcome leakage, but the same underlying donor pools and some regions recur in validation. Monte Carlo intervals do not capture donor-panel uncertainty, tuning uncertainty or the distribution of future genomic regions. The 100-draw frequency grid produces ties, and B200 sensitivity is not independently recalibrated. Secondary Bayesian methods use distinct latent inclusion and coefficient models. Runtime is implementation- and hardware-dependent.

No ordinary-bootstrap stability-selection FDR theorem is claimed. This experiment supplies neither a genome-wide threshold nor a biological-discovery result, and it does not establish a universally superior replacement for Manc-COJO.

## Conclusion

An R bootstrap LASSO selector can replace the stepwise stage in a Manc-COJO workflow while retaining the official joint-estimation implementation. In this experiment it improved recovery of heterogeneous and three-causal signals, but reduced one-causal and shared-effect recovery, cost additional runtime, and did not consistently improve perturbation stability. The reviewed wrapper and paired validation make these tradeoffs reproducible; they do not establish a generally superior replacement or a biological finding.

## Reproducibility and availability

All project files are stored locally on Toshiba. The source, audited results, figures and Word manuscript are available in the public [GitHub repository](https://github.com/Vusumuzi2Mabasa/Multiancestry_Bootstrap_LASSO). The predecessor repository is preserved. The [workflow](workflow.txt) gives commands for provenance verification, numerical controls, independent tuning, threshold freezing, validation, auditing and manuscript generation. Invoke every Rscript with `--vanilla`; personal R startup files are never read.

Public-reference genotype derivatives stay in ignored `data/local/`; raw genotypes and donor identifiers are not committed. [Provenance](results/provenance.json) records exact input hashes and the official executable identity. Reproduction requires those verified local inputs or regeneration using the pinned source study's public-data preparation workflow. [Fit statuses](results/fit_status.tsv), [audit](results/audit.json), [thresholds](results/threshold_summary.tsv) and all per-replicate tables are retained. The Word manuscript is generated from this README with native editable equations and author Vusumuzi Mabasa.

## References

1. Wang X, Wang Y, Yang J, Visscher PM, Wray NR, Yengo L. Multiancestry conditional and joint analysis applied to GWAS summary statistics. Nature Genetics (2026). [doi:10.1038/s41588-026-02799-1](https://doi.org/10.1038/s41588-026-02799-1). [Official source](https://github.com/light156/multi-ancestry-COJO/tree/296953627d2d2b7bc3bfa04cf65524c5013c3b41) and [command interfaces](https://light156.github.io/multi-ancestry-COJO-docs/options/).
2. Bach FR. Bolasso Model consistent Lasso estimation through the bootstrap. ICML (2008). [doi:10.1145/1390156.1390161](https://doi.org/10.1145/1390156.1390161). This study uses thresholded parametric-bootstrap frequencies, not the exact Bolasso support-intersection procedure or a claimed transfer of its asymptotic theorem.
3. Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software (2010). [doi:10.18637/jss.v033.i01](https://doi.org/10.18637/jss.v033.i01). [glmnet documentation](https://glmnet.stanford.edu/articles/glmnet.html).
4. The 1000 Genomes Project Consortium. A global reference for human genetic variation. Nature (2015). [doi:10.1038/nature15393](https://doi.org/10.1038/nature15393). [IGSR data reuse](https://www.internationalgenome.org/IGSR_disclaimer/).
5. Rossen J, Shi H, Strober BJ, et al. MultiSuSiE improves multi-ancestry fine-mapping in All of Us whole-genome sequencing data. Nature Genetics (2026). [doi:10.1038/s41588-025-02450-5](https://doi.org/10.1038/s41588-025-02450-5). Used as an external secondary implementation; no All of Us data are analysed here.
'''
    (ROOT/'README.md').write_text(text,encoding='utf-8',newline='\n')
    print('Research README generated from audited results')
if __name__=='__main__':main()
