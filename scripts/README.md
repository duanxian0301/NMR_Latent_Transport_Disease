# Script provenance and execution order

These files are sanitized copies of the scripts used in the frozen project. Path changes only replace local absolute paths with environment variables; statistical thresholds, seeds, factor sets and analysis logic are unchanged.

Suggested source-chain order:

1. `input_preparation/`: derive the common harmonized SNP panel and per-trait sample-size table after the original public GWAS files have been downloaded.
2. `core/run_exp52_9_gfa.R`: construct the 61-trait matrices, residual-correlation matrix and UK Biobank GFA fit.
3. `core/run_exp52_9_gls_projection.R`: derive resource-specific GLS factor projections and cross-resource summaries.
4. `independent_replication/`: fit Karjalainen GFA independently and compare loading structures.
5. `single_trait_baseline/`: audit and summarize the 61 same-metabolite cross-resource LDSC analyses.
6. `audits/`: export frozen loadings and classify factor associations relative to their constituent metabolites.
7. `locus_lava/`: construct compact locus summaries and run the cross-resource and UKB inter-factor LAVA analyses.

The original download/filter launchers, temporary diagnostics, machine-specific wrappers and scripts that deleted large intermediate downloads after filtering are intentionally excluded. See `PUBLIC_ARCHIVE_QA_CN.md`.

