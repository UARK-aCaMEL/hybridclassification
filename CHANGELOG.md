# aCaMEL/hybridclassification: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### `Fixed`

- The run stopped, without a report, when no sample was classified as a hybrid above `--prob_threshold` ([#8](https://github.com/UARK-aCaMEL/hybridclassification/issues/8)). The report is now built, with empty hybrid and mask lists.
- With `--run_bgc`, a pair with fewer than two candidate hybrids stopped the run in `BGC_HM`. Such pairs now skip genomic clines with a warning.
- `MASK_SAMPLES` failed when none of a hybrid class's calls had a hybrid index (for example, when all of them were parental reference samples).

## v1.0.0 - Schism - [30-09-2026]

Initial release of aCaMEL/hybridclassification, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `--seed` for reproducible runs. It seeds SNP thinning, every ADMIXTURE run, the hybrid simulations, NewHybrids and bgchm. If unset, a seed is derived from the session ID; it is kept by `-resume`, printed at startup and recorded in the report.
- evalAdmix model fit for every K, with a report section for K = 2 (the model used to find candidate hybrids).
- All pipeline parameters in `nextflow_schema.json`, so they are validated and listed by `--help`.
- Documentation: README, usage, output and citations.

### `Fixed`

- DISTRUCT re-ran on every `-resume` (it wrote into CLUMPAK's output directory), and so did everything after it.
- NewHybrids re-ran on every `-resume` (its seeds were drawn at random in the script block).
- `TOP_LOCI_FST` dropped the header line of the per-locus F<sub>ST</sub> table.
- `--geo_data_dir` defaulted to a path relative to the launch directory.
- Only biallelic SNPs are kept. SNPio 1.7 keeps every allele at multiallelic sites, which the NewHybrids, simulation and bgc converters do not handle.
- Strict `nextflow lint` errors.
- The `test_full` profile was the nf-core template's placeholder; it now runs the bundled test data with the default analysis settings.
- The VCF decompressed for AdmixPipe was published to `<pair>/admixpipe/decompress_vcf/`.
- A pair could be analysed with one of its species missing. When the retained loci were genotyped mostly in the more numerous species, the per-sample missing-data filter removed every sample of the other: on the test data, all 56 CHRERY samples in CAMANO_CHRERY, so ADMIXTURE and NewHybrids compared two groups of CAMANO. `SNPIO_FILTER` now applies the per-species missing-data filter (`--pop_cov`) first, and stops with an error if fewer than `--min_species_samples` (default 5) samples of either species remain.

### `Dependencies`

- nf-core template 2.14.1 → 4.1.0
- SNPio 1.3.21 → 1.7.6
- AdmixPipe 3.2 → 3.2.2
- nf-core `tabix/bgzip` and `tabix/tabix` (deprecated) replaced by `htslib/bgziptabix`

### `Deprecated`
