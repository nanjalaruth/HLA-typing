# HLA-typing

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

## Introduction

The pipeline does Human Leukozyte Antigen (HLA) typing using [HLA-VBSEQ]() from high-throughput sequencing data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker and singularity containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility 

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nanjalaruth/HLA-typing -profile test,<docker/singularity/conda/institute>
    ```

4. Start running your own analysis!

    ```bash
    nextflow run nanjalaruth/HLA-typing -profile <docker/singularity/conda/institute> --input '*_R{1,2}.fastq.gz'
    ```
