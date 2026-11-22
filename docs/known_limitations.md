# Known limitations

1. The metabolite panel is enriched for lipid and lipoprotein measures, so latent directions reflect both biology and measurement redundancy.
2. The 13-factor set is conditional on cross-resource QC that used Karjalainen information.
3. Fixed-loading projections and independently estimated factors are different estimands.
4. Same-GWAS comparisons between factors and metabolites are correlated and do not provide an independent test of statistical gain.
5. Most factor-associated lead regions overlap genome-wide significant single-metabolite regions.
6. Global independent-factor correspondence did not exceed the prespecified residual-correlation null, despite close recovery of selected profiles such as the F17 branched-chain amino-acid direction.
7. Inter-factor LAVA is based on UK Biobank factor GWAS only.
8. Original GWAS files, UK Biobank controlled data, LD references, RDS fits and serialized matrix caches are not redistributed.
9. The full primary pipeline therefore requires users to obtain third-party inputs and reconstruct the documented directory structure.
10. NDD, Aging and MVMR feasibility work is outside the positive manuscript result set and is intentionally absent from this archive.

