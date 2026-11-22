# Analysis overview

## 1. Harmonized GWAS panel

Sixty-one NMR metabolite traits were matched between UK Biobank and Karjalainen using GWAS Catalog accessions. Association statistics were aligned by rsID and alleles. The final common panel contained 1,202,568 variants; the bulk variant-level files are not redistributed.

## 2. UK Biobank GFA

GFA was fitted to LD-pruned UK Biobank association Z scores using the regularized LDSC residual-correlation matrix. The frozen primary seed was 101, with seed 202 used for stability checks. Model selection yielded 30 factors. The complete 61 × 30 normalized loading matrix is in `derived/ukb_gfa_loadings_61x30.tsv`.

## 3. Fixed-loading GLS projection

The UK Biobank loading directions were frozen. Genome-wide factor association statistics were computed separately in UK Biobank and Karjalainen using each resource's own regularized residual-correlation matrix. These are transferred projections, not independent factor estimates.

## 4. Frozen factor QC

Thirteen factors were retained for principal cross-resource summaries. The selection used calibration, heritability and directional-concordance criteria and is therefore conditional on cross-resource information. The original 30-factor classifications are preserved in `derived/factor_qc_decision.tsv`.

## 5. Single-trait baseline

Same-metabolite cross-resource LDSC was evaluated for all 61 traits. This provides the reference against which the fixed-loading factor projections were described. It is not an independent statistical test of factor superiority because comparisons share underlying GWAS inputs.

## 6. Independent Karjalainen GFA

Karjalainen GFA was fitted independently and selected 15 factors. One-to-one and nearest-profile loading matches, subspace statistics and null-calibrated comparisons are reported separately from the fixed-loading projection analyses.

## 7. Locus and LAVA analyses

Factor-associated leads were compared with constituent single-metabolite signals. LAVA evaluated same-factor cross-resource local sharing and UKB inter-factor local correlations. Inter-factor LAVA was not replicated in Karjalainen.

## 8. Figures

Figures 2–4 are deterministic transforms of the committed source TSV files. Their plotting scripts do not read the bulk GWAS inputs.

