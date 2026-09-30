# aCaMEL/hybridclassification: Output

## Introduction

This document describes the output produced by the pipeline. The main result for each species pair is an interactive report; the other files are the intermediate and detailed results behind it.

The directories listed below will be created in the results directory after the pipeline has finished. Each species pair in `--combinations` gets its own directory, named `<pair>` (for example `CAMANO_CAMOLI`), and the files inside it are prefixed with the pair name. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Reports](#reports): one interactive report per species pair
- [SNPio](#snpio): subsetting and filtering of SNPs and samples
- [ADMIXTURE](#admixture): ancestry estimation, replicate alignment and model fit
- [Candidate hybrids](#candidate-hybrids): assignment of samples to parental groups and candidates
- [Locus panel](#locus-panel): the most differentiated loci between the parental groups
- [Simulations](#simulations): simulated hybrids and the NewHybrids power analysis
- [NewHybrids](#newhybrids): classification of the samples
- [Triangle plots and masking](#triangle-plots-and-masking): hybrid index, heterozygosity and outlier masking
- [Genomic clines](#genomic-clines): bgchm results, with `--run_bgc`
- [Pipeline information](#pipeline-information): report metrics generated during the workflow execution

### Reports

<details markdown="1">
<summary>Output files</summary>

- `reports/`
  - `<pair>_multiqc_report.html`: the report for each species pair, a standalone HTML file that can be viewed in your web browser.
- `<pair>/report/`
  - `multiqc_report.html`: the same report before the aCaMEL logo is added.
  - `multiqc_data/`: data behind the report tables and plots.
- `<pair>/plot/`, `<pair>/nh/`: the report's figures and tables as separate files.

</details>

Each report, built with [MultiQC](http://multiqc.info), has these sections:

- **ADMIXTURE Results:**
  - the K = 2 ancestry barplot, whose hover text gives each sample's group (P0, P1 or candidate);
  - the evalAdmix fit of K = 2.
- **NewHybrids Results:**
  - posterior class probabilities for each sample, and summary tables by species and site, with and without masking;
  - simulation accuracy, triangle plots, outlier diagnostics and the MCMC trace;
  - with `--site_coords`, a map of the classifications at each site.
- **Genomic Cline Results** (with `--run_bgc`): hybrid indices, genomic clines and the Stan summary.

The report also records the software versions, the value of every parameter and the random seed.

### SNPio

<details markdown="1">
<summary>Output files</summary>

- `<pair>/snpio/snpio_select/`
  - `<pair>.subset.vcf.gz`: the input VCF subset to the samples of the two species, with its `.tbi` index.
- `<pair>/snpio/snpio_filter/`
  - `<pair>.filter.vcf.gz`: the filtered VCF used by ADMIXTURE, NewHybrids and bgchm, with its `.tbi` index.
  - `<pair>_output/`: SNPio's missing-data tables (`reports/missingness/`), filtering Sankey diagram (`plots/`) and logs.
- `<pair>/snpio/snpio_popfilter/`
  - `<pair>.temp.filtered.vcf.gz`: the filtered VCF after a second missing-data filter within the P0, P1 and candidate groups, used for the locus panel and bgchm.

</details>

[SNPio](https://github.com/btmartin721/SNPio) filters SNPs and samples by missing data (within each species, per SNP and per sample), minor allele frequency and physical distance, and keeps only biallelic SNPs. See [Filtering](usage.md#filtering).

### ADMIXTURE

<details markdown="1">
<summary>Output files</summary>

- `<pair>/admixpipe/admixturepipeline/`
  - `<pair>.<K>_<rep>.Q`, `<pair>.<K>_<rep>.P`: ancestry proportions and allele frequencies for each K and replicate.
  - `<pair>.<K>_<rep>.stdout`: ADMIXTURE log, including the cross-validation error and log-likelihood.
  - `<pair>_inds.txt`, `<pair>_pops.txt`: sample and population order of the Q files.
- `<pair>/admixpipe/clumpak/clumpakOutput/`: replicates aligned and grouped into modes for each K, with a barplot of each major mode.
- `<pair>/admixpipe/distruct/`: barplots (`K<K>.pdf`) and the CLUMPP-aligned results of the major mode for each K.
- `<pair>/admixpipe/cvsum/`: cross-validation error and log-likelihood across K.
- `<pair>/admixpipe/evaladmix/`: evalAdmix residual correlations (`.corres`) and heatmaps for each K and replicate.
- `<pair>/admixpipe/bestk/`
  - `<pair>_k2_clumpp_indfile.out`: the aligned K = 2 ancestry proportions used to find candidate hybrids.
  - `<pair>_bestK.txt`, `<pair>_best_clumpp_indfile.out`: the K (from 2 to `--maxk`) with the lowest cross-validation error, and its ancestry proportions.

</details>

[ADMIXTURE](https://dalexander.github.io/admixture/) is run by [AdmixPipe](https://github.com/stevemussmann/admixturePipeline) for K = 1 to `--maxk`, with replicates aligned by [CLUMPAK](https://clumpak.tau.ac.il/). [evalAdmix](https://github.com/GenisGE/evalAdmix) tests how well each K fits: residual correlations near zero mean a good fit, and positive correlations within a group mean the model does not capture its structure.

### Candidate hybrids

<details markdown="1">
<summary>Output files</summary>

- `<pair>/find/`
  - `<pair>_popmap.tsv`: each sample's group from K = 2: `P0`, `P1`, or `ADMIX` (candidate hybrid).

</details>

A sample whose larger K = 2 ancestry proportion is below `--ancestry_threshold` is a candidate hybrid. The others form the parental groups P0 and P1.

### Locus panel

<details markdown="1">
<summary>Output files</summary>

- `<pair>/selected_loci/`
  - `<pair>_fst_per_locus.tsv`: the loci in the panel and their Weir & Cockerham F<sub>ST</sub> between P0 and P1, highest first.
  - `<pair>_top<N>_fst.vcf`: the panel's genotypes, for NewHybrids and the simulations.

</details>

The panel is the `--panel_size` loci with the highest F<sub>ST</sub> between P0 and P1, computed with [VCFtools](https://vcftools.github.io/). If fewer loci pass filtering, all of them are used.

### Simulations

<details markdown="1">
<summary>Output files</summary>

- `<pair>/power_analysis/simulated_data/`
  - `<pair>_simulation.vcf`, `<pair>_simulation.tsv`: simulated individuals of each class, and their classes.
  - `<pair>_NewHybrids.txt`, `<pair>_NH_index_map.tsv`: the NewHybrids input, with real and simulated samples, and the sample behind each NewHybrids index.
- `<pair>/power_analysis/newhybrids/`: NewHybrids results for that input (as in [NewHybrids](#newhybrids)).

</details>

Pure, F1, F2 and backcross individuals are simulated from the P0 and P1 allele frequencies at the panel loci. The share of each simulated class that NewHybrids assigns correctly is shown in the report's _Simulation Assignment Accuracy_ plot. If simulated hybrids of a class are often misassigned, real hybrids of that class cannot be identified reliably with this panel.

### NewHybrids

<details markdown="1">
<summary>Output files</summary>

- `<pair>/newhybrids/`
  - `aa-PofZ.txt`: each sample's posterior probability of each class (P0, P1, F1, F2, Bx0, Bx1).
  - `<pair>_NH_index_map.tsv`: the sample behind each NewHybrids index in `aa-PofZ.txt`.
  - `pi_trace.tsv`: the MCMC trace of the class proportions.
  - `<pair>_NewHybrids.txt`: the NewHybrids input. Samples in P0 and P1 are marked as known parental samples.
  - `aa-*`, `EchoedGtypData.txt`, `<pair>_nh_results.txt`: other NewHybrids output and its log.
- `<pair>/nh/`
  - `<pair>_nh_hybrids.txt`, `<pair>_nh_hybrids_masked.txt`: samples classified as hybrids (F1, F2, Bx0 or Bx1 with posterior above `--prob_threshold`), before and after masking.
  - `<pair>_points.tsv`: counts of each class per site, for the map.

</details>

[NewHybrids](https://github.com/eriqande/newhybrids) assigns each sample a posterior probability of belonging to each of six genotype frequency classes, from the panel loci. A sample is assigned to a class when that probability exceeds `--prob_threshold`.

### Triangle plots and masking

<details markdown="1">
<summary>Output files</summary>

- `<pair>/triangle/`
  - `<pair>_hindex.tsv`: hybrid index and interspecific heterozygosity of every real and simulated sample. It uses the loci whose allele frequencies differ between P0 and P1 by at least `--af_dist_min`; the header gives the number of loci.
  - `<pair>_hindex_fixed.tsv`: the same, from fixed differences only (loci with an allele frequency difference of at least 0.999), if there are any.
  - `<pair>_classification_popmap.tsv`: the group of each real and simulated sample.
- `<pair>/mask/`
  - `<pair>_mask_table.tsv`: for each sample classified as a hybrid, its Mahalanobis distance (D²) and p-value relative to the simulated individuals of its class, and whether it falls inside the ellipse.
  - `<pair>_masked_samples.txt`: the samples outside it, which are masked.

</details>

Each hybrid class occupies a characteristic region of a triangle plot of hybrid index against interspecific heterozygosity ([Fitzpatrick 2012](https://doi.org/10.1186/1471-2148-12-131)). A hybrid call is masked when the sample falls outside the 1 − `--outlier_alpha` ellipse of the simulated individuals of its class. For example, a sample NewHybrids calls F1 but whose heterozygosity is far below that of simulated F1s is masked.

### Genomic clines

<details markdown="1">
<summary>Output files</summary>

With `--run_bgc` only.

- `<pair>/vcf2bgc/`: genotype-likelihood inputs for P0, P1 and the candidate hybrids, with the sample and locus order.
- `<pair>/bgc/`
  - `<pair>_cline_parameter_summary.tsv`: for each locus, the posterior median and credible interval of the cline centre and gradient, and whether its cline is an outlier (e.g. shifted towards P0, or steeper than average).
  - `results_<pair>/text/`: hybrid indices, cline parameters, MCMC draws, sampler diagnostics and Stan summaries.
  - `results_<pair>/plots/`: PDF plots of the hybrid indices, genomic clines and triangle plot.
  - `results_<pair>/rdata/`: the full R results.

</details>

[bgchm](https://github.com/zgompert/bgc-hm) estimates each candidate hybrid's hybrid index, and a genomic cline for each locus that describes how its ancestry changes with the hybrid index. The report's Stan summary gives effective sample sizes and R̂ for checking convergence.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.html`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `hybridclassification_software_mqc_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameters are used when running the pipeline.
  - Parameters used by the pipeline run: `params_<timestamp>.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
