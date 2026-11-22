# Pleiotropic genetic patterns of circulating metabolites

This repository contains the reproducibility code, compact derived results and figure source data for the study **“Pleiotropic genetic patterns of circulating metabolites and their reproducibility across GWAS resources.”**

The study applies genetic factor analysis (GFA) to 61 harmonized nuclear magnetic resonance metabolite GWAS from UK Biobank and evaluates the resulting genetic profiles in the independent Karjalainen GWAS resource. The public archive separates three distinct questions:

1. the composition of UK Biobank GFA loading profiles;
2. cross-resource reproducibility after projecting fixed UK Biobank loadings with resource-specific generalized least squares (GLS); and
3. correspondence with factors estimated independently in Karjalainen.

These are different estimands. A reproducible fixed-loading projection does not establish recovery of the same independently estimated latent factor.

## Study scope

- 61 harmonized NMR metabolites;
- 30 UK Biobank GFA factors;
- 13 factors retained under the frozen cross-resource QC rule;
- 15 factors estimated by independent Karjalainen GFA;
- resource-specific GLS projections;
- single-metabolite cross-resource LDSC reference analyses;
- compact locus-level and LAVA summaries;
- source data and plotting code for manuscript Figures 2–4.

The archive does **not** present NDD, Aging or MVMR feasibility analyses as positive manuscript results.

## Repository map

| Path | Contents |
|---|---|
| `scripts/` | Sanitized copies of the scripts used for input preparation, GFA, GLS projection, independent-factor comparison, single-trait baselines, locus summaries and LAVA |
| `metadata/` | The 61-trait crosswalk, public accessions, input manifest, column dictionaries and provenance |
| `derived/` | Small final loading matrices and compact result summaries; no bulk GWAS files |
| `figures/` | Exact Figure 2–4 source TSV files and plotting scripts |
| `validation/` | Release checks, expected hashes and manuscript/supplement numeric QA |
| `environment/` | Python and R software requirements |
| `config/` | Example path configuration for external resources |
| `docs/` | Analysis overview, decision definitions, limitations and availability notes |

## Data sources

The original metabolite association statistics are not redistributed. Public GWAS Catalog accessions for both resources are listed in [`metadata/trait_crosswalk_61.tsv`](metadata/trait_crosswalk_61.tsv) and as stable links in [`metadata/input_manifest.tsv`](metadata/input_manifest.tsv). Readers should obtain the original summary statistics from the named repositories and comply with the source terms.

Individual-level UK Biobank data, controlled resources, bulk 1.2-million-variant metabolite or factor GWAS files, LD reference panels and locally cached objects are excluded.

## Reproduce Figures 2–4

Create the Python environment described in [`environment/environment.yml`](environment/environment.yml), then run from the repository root:

```bash
python figures/make_figure2.py
python figures/make_figure3.py
python figures/make_figure4.py
```

Each script reads the adjacent source TSV and writes PDF, SVG and 600-dpi PNG files. The committed source-data hashes are checked by:

```bash
python validation/check_release.py
```

## Re-running analytical stages

The analytical scripts preserve the frozen seeds, thresholds and selection rules used in the study. They require external GWAS summary statistics and reference resources that are not bundled here. Configure their locations with the environment variables documented in [`config/config.example.yaml`](config/config.example.yaml). The scripts are published primarily for provenance and method transparency; they are not a promise that restricted or third-party inputs can be redistributed.

Key software includes R 4.4.0, GFA 1.0.0.0470, PLINK 1.9, LDSC 1.0.1 and LAVA 0.1.5. See [`environment/software_versions.md`](environment/software_versions.md).

## Interpretation boundaries

- The 13-factor analysis is conditional on the frozen cross-resource QC selection.
- Factor orientation is arbitrary; sign alone is not a biological mechanism.
- Most factor–lead associations overlap genome-wide significant single-metabolite regions. “Multivariate-threshold candidate” does not mean an LD-independent novel locus.
- Inter-factor LAVA was conducted in UK Biobank only.
- The independent Karjalainen factor space did not show global excess correspondence under the prespecified null, although specific profiles such as the F17 branched-chain amino-acid direction were closely recovered.

See [`docs/known_limitations.md`](docs/known_limitations.md) for the full boundaries.

## Citation

The manuscript author list, affiliations and preferred repository citation are still pending author completion. Until these fields are verified, cite the repository URL and the manuscript title above. A `CITATION.cff` and Zenodo metadata record will be added only after the authors approve the final metadata.

## Licensing

Source code is released under the [MIT License](LICENSE). Original documentation, metadata, figure source tables and derived summary tables are released under [CC BY 4.0](LICENSE-DATA). Third-party GWAS, UK Biobank controlled resources, LD reference data and other source datasets are not redistributed and are not covered by these repository licences.

