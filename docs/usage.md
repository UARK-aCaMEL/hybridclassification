# aCaMEL/hybridclassification: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files. Run `nextflow run UARK-aCaMEL/hybridclassification --help` to list them._

## Introduction

aCaMEL/hybridclassification detects and classifies hybrids between pairs of species. It needs a VCF, a species map, a population map and a list of species pairs. Site coordinates and map layers are optional.

Each species pair is analysed separately. The pipeline subsets the VCF to the two species and filters it, then runs ADMIXTURE and splits the samples:

- **Candidate hybrids:** samples with mixed ancestry at K = 2.
- **Parental reference groups:** all the others, labelled P0 and P1.

Next it picks the loci most differentiated between P0 and P1. It uses them to simulate hybrids of known class and to classify the candidates with NewHybrids.

## Inputs

All the text inputs are tab-delimited, with no header. Sample IDs must match the VCF header.

### VCF

A single multi-sample VCF of SNPs, as `.vcf` or `.vcf.gz`. The pipeline indexes it, and keeps only biallelic SNPs.

```bash
--input '[path to VCF file]'
```

### Species map

The first column is the sample ID. The second is the species.

```bash
--speciesmap '[path to species map]'
```

```text title="speciesmap.tsv"
Z01CAMANO01	CAMANO
Z01CHRERY01	CHRERY
Z03CAMOLI01	CAMOLI
```

### Population map

The first column is the sample ID. The second is the population or sampling site. Sites are used to summarise the classifications and to place samples on the map.

```bash
--popmap '[path to popmap file]'
```

```text title="popmap.tsv"
Z01CAMANO01	Z01
Z01CHRERY01	Z01
Z03CAMOLI01	Z03
```

### Species pairs

Two species IDs from the species map per line. Each line is analysed separately, with its own subdirectory and report.

```bash
--combinations '[path to pairs file]'
```

```text title="combinations.tsv"
CAMANO	CAMOLI
CAMANO	CHRERY
```

### Site coordinates (optional)

Site ID (as in the population map), latitude and longitude, in decimal degrees. When provided, the report includes a map of the NewHybrids classifications at each site.

```bash
--site_coords '[path to coordinates file]'
```

```text title="site_coords.tsv"
Z01	35.90463	-91.63537
Z02	36.09227	-91.75397
```

### Map layers (optional)

Extra vector layers, such as streams or range boundaries, can be drawn on the map. Put the layer files in a directory, and describe them in a JSON file:

```bash
--geo_data_dir '[directory of layer files]' --geo_data_config '[layers JSON]'
```

```json title="layers.json"
[
  {
    "path": "test_geo_data/streams.shp",
    "z_order": 1,
    "style": { "color": "#0066ff", "weight": 1, "opacity": 1.0 }
  }
]
```

Each `path` starts with the name of the `--geo_data_dir` directory. `style` takes [Leaflet path options](https://leafletjs.com/reference.html#path-option). See [`assets/test_geo_data.json`](../assets/test_geo_data.json) for an example.

## Main settings

### Filtering

For each pair, SNPio applies these filters in order:

1. Per-species missing data (`--pop_cov`).
2. Monomorphic and multiallelic sites removed.
3. Per-SNP missing data (`--snp_cov`).
4. Minor allele frequency (`--min_maf`).
5. Per-sample missing data (`--ind_cov`).
6. Thinning to one SNP per `--thin_dist` bp. The SNP kept in each window is chosen at random, from `--seed`.

The per-species filter comes first, so the SNPs kept are genotyped in both species. If fewer than `--min_species_samples` (default 5) samples of either species remain, the pipeline stops with an error that names the species. Otherwise ADMIXTURE would split the one remaining species into two groups, and NewHybrids would call "hybrids" between them. To continue, relax the missing-data thresholds or drop the pair from `--combinations`.

### Candidate hybrids and parental groups

ADMIXTURE runs for K = 1 to `--maxk` (default 2), with replicates aligned by CLUMPAK. Candidate hybrids come from K = 2. A sample whose larger ancestry proportion is below `--ancestry_threshold` (default 0.9) is a candidate hybrid; the rest are assigned to P0 or P1 by their main cluster.

P0 and P1 are used as the known parental references in NewHybrids, so check that the K = 2 clusters really correspond to the two species. The report's ADMIXTURE barplot and evalAdmix residuals help with this.

Advanced-generation hybrids that breed among themselves can form their own cluster in ADMIXTURE. Such a cluster would then be treated as a parental reference.

### NewHybrids and the simulations

The loci passed to NewHybrids are the `--panel_size` (default 500) with the highest Weir & Cockerham F<sub>ST</sub> between P0 and P1. Loci that do not vary in the pair (F<sub>ST</sub> undefined) come last.

**Simulations:** from the P0 and P1 allele frequencies at those loci, the pipeline simulates `--sample_size` (default 10) individuals of each class, in `--n_reps` (default 4) replicates. The classes are pure P0, pure P1, F1, F2, backcross to P0 and backcross to P1. NewHybrids is run on the simulated and real samples together, and the report shows how often each simulated class is assigned correctly. This is the pipeline's measure of power.

**Classification:** in the real data, a sample is assigned to a class when its NewHybrids posterior probability exceeds `--prob_threshold` (default 0.9); otherwise it is unassigned. `--nh_burnin` and `--nh_sweeps` set the length of the NewHybrids MCMC run. Check the trace plot in the report.

The panel is chosen from the same parental samples that the simulations are drawn from. That can make power look higher than it would be for independent data (high-grading bias).

### Outlier masking

The pipeline also computes, for each sample, the hybrid index and interspecific heterozygosity (the triangle plot), from the loci whose allele frequencies differ between P0 and P1 by at least `--af_dist_min` (default 0.7).

A hybrid call is masked (treated as unassigned) when the sample falls outside the 1 − `--outlier_alpha` (default 95%) Mahalanobis ellipse of the simulated individuals of its class. The report shows the summary table both with and without masking.

### Genomic clines (optional)

`--run_bgc` estimates hybrid indices and genomic clines with bgchm, for the candidate hybrids against P0 and P1. It uses the loci that pass `--af_dist_min`. `--bgc_iters`, `--bgc_burnin` (the warm-up proportion) and `--bgc_thin` set the MCMC run.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run UARK-aCaMEL/hybridclassification \
    --input genotypes.vcf.gz \
    --speciesmap speciesmap.tsv \
    --popmap popmap.tsv \
    --combinations combinations.tsv \
    --outdir <OUTDIR> \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run UARK-aCaMEL/hybridclassification -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: 'genotypes.vcf.gz'
speciesmap: 'speciesmap.tsv'
popmap: 'popmap.tsv'
combinations: 'combinations.tsv'
outdir: './results/'
<...>
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull UARK-aCaMEL/hybridclassification
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [aCaMEL/hybridclassification releases page](https://github.com/UARK-aCaMEL/hybridclassification/releases) and find the latest pipeline version - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

The version number and the value of every parameter are recorded in the _Workflow Summary_ section of each report. The parameters are also saved to `pipeline_info/params_<timestamp>.json`, which can be passed back with `-params-file` to repeat a run.

> [!TIP]
> If you wish to share such a parameter file (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

The random steps are seeded from `--seed`: SNP thinning, every ADMIXTURE K and replicate, the hybrid simulations, NewHybrids and bgchm. A rerun with the same seed, inputs and parameters repeats them. If `--seed` is not set, a seed is derived from the run's session ID. It differs between runs, is kept when you `-resume`, and is printed at startup and recorded in the report. Rerun with `--seed <value>` to reproduce that run.

Multithreaded ADMIXTURE is not bit-for-bit deterministic, so the K = 2 ancestry proportions, and occasionally the candidate hybrids, can differ slightly between reruns. For exact reproduction, give `ADMIXTUREPIPELINE` one CPU (see [Resource requests](#resource-requests)).

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer) - see below.

> [!IMPORTANT]
> This pipeline requires a container engine. Conda is not supported, because several steps use purpose-built containers.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

- `test`
  - A profile with a complete configuration for automated testing
  - Uses the bundled test data, with short NewHybrids runs and a 100-locus panel, so needs no other parameters
- `test_full`
  - The same test data, run with the default analysis settings
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `emulate_amd64`
  - Use together with `docker` on Apple Silicon and other ARM machines, to run the x86-64 containers under emulation
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)

On the University of Arkansas AHPCC Pinnacle cluster, add the bundled Slurm configuration with `-c ahpcc.config -profile singularity`.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code indicating it ran out of resources, it will automatically be resubmitted once with double the requests. If it fails again, the pipeline execution is stopped.

ADMIXTURE (`ADMIXTUREPIPELINE`) and NewHybrids (`RUN_NEWHYBRIDS`, `POWER_ANALYSIS`) are the most demanding steps. ADMIXTURE's run time grows with `--maxk` and the number of replicates. NewHybrids' run time grows with `--nh_sweeps`, `--panel_size` and, for the power analysis, `--n_reps` × `--sample_size`.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

For example, ADMIXTURE runs 10 replicates per K with 10-fold cross-validation by default (2 and 2 in the `test` profile). This is set by `ext.args` for `ADMIXTUREPIPELINE` in [`conf/modules.config`](../conf/modules.config).

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time. Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
