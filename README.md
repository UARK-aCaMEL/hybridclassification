<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/acamel_hybridclassification_logo_dark.png">
    <img alt="aCaMEL/hybridclassification" src="docs/images/acamel_hybridclassification_logo_light.png" width="400">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/UARK-aCaMEL/hybridclassification/actions/workflows/nf-test.yml/badge.svg)](https://github.com/UARK-aCaMEL/hybridclassification/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/UARK-aCaMEL/hybridclassification/actions/workflows/linting.yml/badge.svg)](https://github.com/UARK-aCaMEL/hybridclassification/actions/workflows/linting.yml)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/UARK-aCaMEL/hybridclassification)

## Introduction

**aCaMEL/hybridclassification** is a bioinformatics pipeline that detects and classifies hybrids between pairs of species from SNP data. It takes a multi-sample VCF, a map of samples to species, and a list of species pairs.

- **Candidate hybrids:** for each pair, it finds candidates with ADMIXTURE.
- **Classification:** it classifies samples as pure, F1, F2 or backcross with NewHybrids.
- **Checks:** it tests how well those classes can be told apart, using simulated hybrids.
- **Genomic clines (optional):** it can also estimate hybrid indices and genomic clines with bgchm.

Each pair gets an interactive report.

1. Subset the VCF to the two species and filter SNPs and samples ([`SNPio`](https://github.com/btmartin721/SNPio))
2. Estimate ancestry with ADMIXTURE, with replicates aligned by CLUMPAK, and assess the fit of K = 2 ([`AdmixPipe`](https://github.com/stevemussmann/admixturePipeline), [`ADMIXTURE`](https://dalexander.github.io/admixture/), [`CLUMPAK`](https://clumpak.tau.ac.il/), [`evalAdmix`](https://github.com/GenisGE/evalAdmix))
3. Flag samples with mixed K = 2 ancestry as candidate hybrids; the others define the two parental groups
4. Choose the loci most differentiated between the parental groups ([`VCFtools`](https://vcftools.github.io/), Weir & Cockerham F<sub>ST</sub>)
5. Simulate pure, F1, F2 and backcross individuals, and measure how often NewHybrids assigns them correctly
6. Classify samples with [`NewHybrids`](https://github.com/eriqande/newhybrids), and flag hybrid calls that do not fit the simulated hybrid index and heterozygosity of their class
7. Optionally, estimate hybrid indices and genomic clines ([`bgchm`](https://github.com/zgompert/bgc-hm))
8. Build an interactive report per species pair ([`MultiQC`](http://multiqc.info/))

## Getting started

You need [Nextflow](https://www.nextflow.io/docs/latest/install.html) (≥ 25.10.4) and a container engine: Docker, Singularity/Apptainer or Podman. Conda is not supported.

To check that everything works, and to see what the outputs look like, run the pipeline on the bundled test dataset. It contains RAD-seq SNPs for 1,573 samples of 15 stream-fish species from 69 sites, and the test runs two species pairs:

```bash
nextflow run UARK-aCaMEL/hybridclassification -profile test,docker --outdir test_results
```

Use `-profile test,singularity` on HPC systems, or add `emulate_amd64` on Apple Silicon (`-profile test,docker,emulate_amd64`). When the run finishes, open the reports in `test_results/reports/` in a web browser.

Nextflow downloads the pipeline from GitHub the first time you run `UARK-aCaMEL/hybridclassification`. To run from a local copy instead, for example to modify the pipeline, clone the repository and run `main.nf`:

```bash
git clone https://github.com/UARK-aCaMEL/hybridclassification.git
cd hybridclassification
nextflow run main.nf -profile test,docker --outdir test_results
```

## Usage

The pipeline needs four input files, all tab-delimited with no header except the VCF:

| Parameter        | Contents                                                                |
| ---------------- | ----------------------------------------------------------------------- |
| `--input`        | A multi-sample VCF (`.vcf` or `.vcf.gz`)                                |
| `--speciesmap`   | Sample ID and species, one sample per line                              |
| `--popmap`       | Sample ID and population (sampling site), one sample per line           |
| `--combinations` | Two species IDs per line; each line is a pair to test for hybridization |

`combinations.tsv`:

```text
CAMANO	CAMOLI
CAMANO	CHRERY
```

Now, you can run the pipeline using:

```bash
nextflow run UARK-aCaMEL/hybridclassification \
   -profile <docker/singularity/.../institute> \
   --input genotypes.vcf.gz \
   --speciesmap speciesmap.tsv \
   --popmap popmap.tsv \
   --combinations combinations.tsv \
   --outdir <OUTDIR>
```

Add `--site_coords` to map the classifications by site, and `--run_bgc` for genomic clines. See the [usage documentation](docs/usage.md) for the input formats and the main settings. Run with `--help` for every parameter.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Pipeline output

The main output is one interactive report per species pair, `<OUTDIR>/reports/<pair>_multiqc_report.html`. For the other output files, see the [output documentation](docs/output.md).

## Credits

aCaMEL/hybridclassification was originally written by [Tyler K. Chafin](https://github.com/tkchafin).

We thank the following people for their assistance in the development of this pipeline:

- [Steven M. Mussmann](https://github.com/stevemussmann), author of AdmixPipe
- [Bradley T. Martin](https://github.com/btmartin721), author of SNPio

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md). Bugs and questions can be reported on the [issue tracker](https://github.com/UARK-aCaMEL/hybridclassification/issues).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use aCaMEL/hybridclassification for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

Please also cite the tools used by the pipeline:

- **SNPio**: Martin BT, Monaco DR, Sharabi N, Mussmann SM, Chafin TK (2026). SNPio: a Python interface for population genomic data processing. _BMC Bioinformatics_. doi: [10.1186/s12859-026-06546-5](https://doi.org/10.1186/s12859-026-06546-5)
- **ADMIXTURE**: Alexander DH, Novembre J, Lange K (2009). Fast model-based estimation of ancestry in unrelated individuals. _Genome Research_ 19:1655–1664. doi: [10.1101/gr.094052.109](https://doi.org/10.1101/gr.094052.109)
- **AdmixPipe**: Mussmann SM, Douglas MR, Chafin TK, Douglas ME (2020). AdmixPipe: population analyses in Admixture for non-model organisms. _BMC Bioinformatics_ 21:337. doi: [10.1186/s12859-020-03701-4](https://doi.org/10.1186/s12859-020-03701-4)
- **AdmixPipe v3**: Mussmann SM, Douglas MR, Chafin TK, Douglas ME (2023). AdmixPipe v3: facilitating population structure delimitation from SNP data. _Bioinformatics Advances_ 3(1):vbad168. doi: [10.1093/bioadv/vbad168](https://doi.org/10.1093/bioadv/vbad168)
- **CLUMPAK**: Kopelman NM, Mayzel J, Jakobsson M, Rosenberg NA, Mayrose I (2015). Clumpak: a program for identifying clustering modes and packaging population structure inferences across K. _Molecular Ecology Resources_ 15:1179–1191. doi: [10.1111/1755-0998.12387](https://doi.org/10.1111/1755-0998.12387)
- **evalAdmix**: Garcia-Erill G, Albrechtsen A (2020). Evaluation of model fit of inferred admixture proportions. _Molecular Ecology Resources_ 20:936–949. doi: [10.1111/1755-0998.13171](https://doi.org/10.1111/1755-0998.13171)
- **NewHybrids**: Anderson EC, Thompson EA (2002). A model-based method for identifying species hybrids using multilocus genetic data. _Genetics_ 160:1217–1229. doi: [10.1093/genetics/160.3.1217](https://doi.org/10.1093/genetics/160.3.1217)
- **hybriddetective** (simulation-based power analysis): Wringe BF, Stanley RRE, Jeffery NW, Anderson EC, Bradbury IR (2017). hybriddetective: a workflow and package to facilitate the detection of hybridization using genomic data in R. _Molecular Ecology Resources_ 17(6). doi: [10.1111/1755-0998.12704](https://doi.org/10.1111/1755-0998.12704)
- **Triangle plots**: Fitzpatrick BM (2012). Estimating ancestry and heterozygosity of hybrids using molecular markers. _BMC Evolutionary Biology_ 12:131. doi: [10.1186/1471-2148-12-131](https://doi.org/10.1186/1471-2148-12-131)
- **bgchm** (with `--run_bgc`): Gompert Z, DeRaad DA, Buerkle CA (2024). A next generation of hierarchical Bayesian analyses of hybrid zones enables model-based quantification of variation in introgression in R. _Ecology and Evolution_ 14(11):e70548. doi: [10.1002/ece3.70548](https://doi.org/10.1002/ece3.70548)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
