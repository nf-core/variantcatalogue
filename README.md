# ![nf-core/variantcatalogue](docs/images/nf-core-variantcatalogue_logo_light.png#gh-light-mode-only) ![nf-core/variantcatalogue](docs/images/nf-core-variantcatalogue_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/variantcatalogue/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/variantcatalogue)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23variantcatalogue-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/variantcatalogue)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/variantcatalogue** is a bioinformatics best-practice analysis pipeline for Workflow designed to generate variant catalogues, a list of variants and their frequencies in a population, from whole genome sequences.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

A. Mapping subworkflow
1. Index the reference genome ([`BWA index`](https://bio-bwa.sourceforge.net/bwa.shtml))
2. Trim the reads ([`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic))
3. Align the reads ([`BWA mem`](https://bio-bwa.sourceforge.net/bwa.shtml))
4. Index the aligned bam files ([`Samtools index`](http://www.htslib.org))
5. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
6. Coverage calculation ([`Mosdepth`](https://github.com/brentp/mosdepth))
7. Alignemnt QC ([`Picard Collect WGS metrics`](https://broadinstitute.github.io/picard/))
8. Alignemnt QC ([`Picard Collect Alignment Summary Metrics`](https://broadinstitute.github.io/picard/))
9. Alignemnt QC ([`Quality Scores distribution`](https://broadinstitute.github.io/picard/))
10. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

B. SNV/Indel subworkflow
1. Variant calling ([`DeepVariant`](https://github.com/google/deepvariant))
2. Joint calling ([`GLnexus`](https://github.com/dnanexus-rnd/GLnexus))
3. Split multi-allelic variants ([`Bedtools norm`](https://bedtools.readthedocs.io/en/latest/))
4. Redefine varaints ID ([`Bedtools Annotate`](https://bedtools.readthedocs.io/en/latest/))
5. Sample QC and sex inference ([`Hail`](https://hail.is))
6. Variant QC and frequency calculation ([`Hail`](https://hail.is))
7. Annotation ([`VEP`](https://useast.ensembl.org/info/docs/tools/vep/index.html))

C. Mitochondrial Subworkflow
1. Index the mitochondrial reference genome ([`BWA index`](https://bio-bwa.sourceforge.net/bwa.shtml))
2. Index the shifted mitochondrial reference genome ([`BWA index`](https://bio-bwa.sourceforge.net/bwa.shtml))
3. Extract the reads mapping to the MT chromosome ([`GATK`](https://gatk.broadinstitute.org))
4. Revert the sam to a fastq file ([`GATK`](https://gatk.broadinstitute.org))
5. Align the reads against the mitochondrial reference genome ([`BWA mem`](https://bio-bwa.sourceforge.net/bwa.shtml))
6. Align the reads against the shifted mitochondrial reference genome ([`BWA mem`](https://bio-bwa.sourceforge.net/bwa.shtml))
7. Index the aligned bam files ([`Samtools index`](http://www.htslib.org))
8. Mark duplicates ([`GATK`](https://gatk.broadinstitute.org))
9. Colect HS metrics ([`Picard`](https://broadinstitute.github.io/picard/))
10. Call variants with Mutect2 ([`GATK`](https://gatk.broadinstitute.org))
11. Liftover the variants mapped against the shifted reference genome ([`GATK`](https://gatk.broadinstitute.org))
12. Merge the VCFs files ([`GATK`](https://gatk.broadinstitute.org))
13. Variant QC and frequency calculation ([`Hail`](https://hail.is))


D. Structural varaint subworkflow
May be included soon


## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a simple three lines command:

   ```
   wget https://raw.githubusercontent.com/scorreard/Variant_catalogue_pipeline/main/testdata/reference/hg38_full_analysis_set_plus_decoy_hla_chr20_X_Y_MT.fa.gz
   gzip -d hg38_full_analysis_set_plus_decoy_hla_chr20_X_Y_MT.fa.gz > hg38_full_analysis_set_plus_decoy_hla_chr20_X_Y_MT.fa
   bash nextflow run nf-core/variantcatalogue -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run nf-core/variantcatalogue --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Documentation

More details about the pipeline can be found in the following preperint [Correard et al., 2022](https://www.biorxiv.org/content/10.1101/2022.10.03.508010v2)

A bytesize talk about the varaint catalogue is available on the [nf-core website](https://nf-co.re/events/2023/bytesize_variantcatalog)

The nf-core/variantcatalogue pipeline comes with documentation about the pipeline [usage](https://nf-co.re/variantcatalogue/usage), [parameters](https://nf-co.re/variantcatalogue/parameters) and [output](https://nf-co.re/variantcatalogue/output).

## Credits

nf-core/variantcatalogue was originally written by @scorreard.

We thank the following people for their extensive assistance in the development of this pipeline:
@WyethWasserman
@melsiddieg
@brittanyhewitson

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#variantcatalogue` channel](https://nfcore.slack.com/channels/variantcatalogue) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use  nf-core/variantcatalogue for your analysis, please cite it using the following doi: [10.1101/2022.10.03.508010](https://doi.org/10.1101/2022.10.03.508010)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
