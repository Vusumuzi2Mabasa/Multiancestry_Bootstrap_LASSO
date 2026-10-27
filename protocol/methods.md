# Selection replacement specification

Author: Vusumuzi Mabasa. All phenotypes are simulated. This is an experimental wrapper around Manc-COJO v1.1.0, not software endorsed by its authors. The original executable is unchanged.

## Retained interfaces and changes

The upstream source at commit 296953627d2d2b7bc3bfa04cf65524c5013c3b41 separates selection (`src/slct_loop.cpp`, `MACOJO::slct_loop`) from the joint-only dispatcher (`src/macojo.cpp`, `MACOJO::entry_function`) and final estimator (`Cohort::calc_joint_effects`). Original selection adds conditional inverse-variance meta-analysis signals and can remove signals through backward selection. Our R selector replaces that entire search. Both pipelines use the original `--effect-mode GCTA` joint estimator, official collinearity checks (default 0.9), allele orientation, and the identical `.ma`, `.bim`, signed `.ld` and `.frq` files. Both apply the same explicit common eligibility mask with `--extract`; reference frequency differences must be strictly below 0.2 in both populations. LASSO passes its proposed set to `--cojo-joint --extract`. There is no subsequent stepwise P-value pruning in the replacement. Official joint estimation can remove collinear variants, which are counted and reported separately.

All eight marginal rows remain in each `.ma` file even when an extract mask is used. This matters: upstream `read_sumstat` computes phenotype variance from all parsed marginal rows, including rows outside the selected set. Truncating the `.ma` files to selected markers would change the joint estimator's scaling. Numerical controls compare the summary reconstruction with official coefficients and SEs, and replay stepwise selected sets through joint-only mode.

## Working summary model

For population k and marker j, let b, s, f and N be the harmonized marginal coefficient, SE, effect-allele frequency and sample size. Reproduce the default GCTA scaling in upstream `src/read_files.cpp`:

$$v_{jk}=2f_{jk}(1-f_{jk}),\quad V_k=\operatorname{median}_j\{v_{jk}N_{jk}[s_{jk}^2+b_{jk}^2/(N_{jk}-1)]\},$$
$$n^{eff}_{jk}=(V_k/v_{jk}-b_{jk}^2)/s_{jk}^2+1,\quad d_{jk}=n^{eff}_{jk}v_{jk}.$$

On the eligible set, Q has diagonal d and off-diagonal entries

$$Q_{ij,k}=R_{ij,k}\min(n^{eff}_{ik},n^{eff}_{jk})\sqrt{v_{ik}v_{jk}},\qquad q_k=\operatorname{diag}(d_k)b_k.$$

The working likelihood is q conditional on beta distributed as normal with mean Q beta and covariance V Q. This covariance is an explicit approximation when Q is reconstructed with external LD, the total-phenotype proxy V, and effective sample sizes. It is not the exact finite-sample distribution of all reported GWAS statistics. No external-reference uncertainty, overlapping-cohort covariance or admixture is integrated. Population cohorts are independent; variants within populations remain correlated. We reject non-positive-definite matrices rather than conceal a repair.

Set H=Q/V, S=diag(sqrt(diag(H))), A=S inverse H S inverse, z=S inverse q/V=S b, and t=S beta. Minimize

$$\sum_{k=1}^2\left\{\frac12t_k^TA_kt_k-z_k^Tt_k+\lambda\sum_j|t_{jk}|\right\}.$$

This is ordinary elementwise LASSO with population-specific coefficients and one shared tuning constant. It separates across populations; it neither pools effect estimates nor imposes shared support. A marker is selected in a replicate if either population has a nonzero coefficient. The allele-scale penalty is lambda times sqrt(d/V) times absolute beta; relative to the average per-population least-squares loss the coefficient threshold is proportional to sqrt(V/d), accommodating unequal information. It is not equal penalty strength per allele and not a group LASSO. Compared with Manc's fixed-effect inverse-variance meta-analysis selection criterion, union selection can retain population-specific or opposite-sign associations. That estimand change is intrinsic to this specified replacement, not a hidden claim of equivalence.

## Parametric bootstrap and lambda selection

For each population and bootstrap draw, generate z star = z + L e, where L L transpose = A and e is a fresh vector of independent standard normals. Equivalently q star is normal with mean q and covariance V Q. The center corresponds to the saturated pilot Q inverse q, which is feasible in these eight-marker positive-definite panels. This perturbs the observed fitted mean rather than generating a global-null bootstrap for every dataset. It can carry noise or external-LD artifacts into the resamples. The bootstrap leaves LD, sample sizes, allele frequencies and V fixed; it is not a participant bootstrap and cannot be applied to dense singular regions without a new specified pilot and regularization analysis.

Fit each bootstrap at the fixed lambda; do not retune lambda inside resamples. Selection frequency is the proportion of 100 union-support selections. Lambda is chosen globally from 1, 1.5, 2 and 2.5 using separate tuning simulations, according to the exact rule in design.json. This tuning uses simulated truth only in the development experiment; evaluation truth never determines lambda or thresholds. Independent null calibration then chooses a stratum-specific frequency cutoff. Given C calibration maxima, use rank ceiling((C+1)0.95) and strict exceedance. Ties can make the procedure conservative; a cutoff of one yields an empty selector and is not overridden. The final stage uses C=199 and independent evaluation streams. The same rank rule calibrates original Manc first-entry scores and secondary Bayesian maximum PIPs. This is empirical locus global-null calibration, not nonnull FDR control or a genome-wide significance recommendation.

## Diagnostics and evaluation

R vectorized coordinate descent solves many bootstrap columns simultaneously; every fit must pass the KKT residual tolerance of 1e-6. Controls compare glmnet on Cholesky pseudo-data, an orthogonal soft-threshold solution, allele flips, marker permutations, and Gaussian covariance. Separate genuine within-population row-bootstrap and fixed-design Gaussian residual simulations use regenerated participant data; these are diagnostics, not the primary method. Row resampling changes the empirical design while the primary parametric bootstrap holds it fixed, so identical frequencies are not expected.

All methods are evaluated on identical fresh simulated datasets, with both matched external LD and an in-sample LD diagnostic. Six fixed regions and CEU paired with YRI or LWK give 12 strata. Each stratum has 60 evaluation nulls and 30 of each of four alternatives, separate from 199 calibration nulls. Genotypes pair whole local haplotypes from the original disjoint study-generating pools; all phenotypes are newly simulated with noise SD one. The 2000 versus 250 sample-size imbalance and original one-, shared-, heterogeneous- and three-causal architectures are retained.

Primary coefficient comparisons use official joint estimates on the respective final selected sets, assigning zero elsewhere. These refits are descriptive post-selection estimates: official SEs are checked numerically but are not reported as valid post-selection confidence intervals. Secondary exact and MultiSuSiE models use their posterior mean estimators; selection scores are masked for ineligible variants, but their latent model retains all eight markers. This secondary-panel difference is disclosed rather than conflated with the primary controlled comparison.

Exact and one-to-one LD-aware recall use the original causal union as denominator, including causes removed by filtering. LD-aware matching uses study r squared at least 0.8 and is reported for each population. Error metrics include null any-selection probability, nonnull mean FDP, per-population all-marker MSE, and bias of the African-source effect 0.28. Monte Carlo intervals use fixed-stratum averages and within-dataset paired differences. They condition on the frozen calibration thresholds and donor panels.

Algorithmic repeat stability compares independent B100 draws on unchanged data. Common perturbation stability compares both primary pipelines after adding one identical working-model perturbation to marginal coefficients, holding marginal SEs and frequencies fixed and reconstructing each pipeline's scaling again. This is a defined summary perturbation diagnostic, not a new independent cohort or a calibrated bootstrap confidence region. Empty-empty Jaccard is one; null and signal strata are reported separately. B200 stability averages the two frequency runs and uses the same B100 cutoff only as a numerical sensitivity diagnostic, without claiming B200 error calibration.

## Scope and negative findings

No stability-selection FDR guarantee, Bayesian PIP interpretation or confidence-interval validity attaches to bootstrap frequencies. Small donor pools, sparse eight-marker regions and reuse of known development regions limit generalization. External LD error can affect both the selection objective and retained joint estimator. Failure or nonconvergence is recorded, never silently dropped; numerical or calibration failures stop the shard. New production results require the protocol, implementation and frozen choices to agree by hash. Results that favor the original method or reveal weak recovery must remain in the report.
