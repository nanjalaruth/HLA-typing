# HLA-typing

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

## Introduction

The pipeline does Human Leukocyte Antigen (HLA) typing using [HLA-VBSEQ](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-16-S2-S7) from high-throughput sequencing data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker and singularity containers making installation trivial and results highly reproducible.

## Installation 
1. [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
2. [Docker](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04) 
3. [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)

## Running the pipeline
The pipeline does not require installation as `NextFlow` will automatically fetch it from `GitHub`.

### Test data
To execute the pipeline on test dataset run:

    ```bash
    nextflow run nanjalaruth/HLA-typing -profile test,slurm,<docker/singularity/conda/institute>
    ```
### Own data
Start running your own analysis either by using flags as shown below:

    ```bash
    nextflow run main.nf -profile singularity -resume --input "*_{1,2}*.fastq.gz" \
    --reference_genome "path to the human reference genome <hg19>" \
    --hla_ref "path to the hla reference genome <hla_all_v2.fasta>" \
    --hla_txt_file "path to the <Allelelist_v2.txt> file"
    ```
 or run your own analysis by modifying the test.config file to suit the path to your data location and then run the command as below:
 
    ```
    nextflow run main.nf -profile singularity -c <path to your edited config file> -resume
    ```

## Arguments

### Required Arguments
| Argument  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| -profile  | \<base,slurm\>                    | Configuration profile to use.                                       |
| --input  | \</project/\*\_{R1,R2}\*.fastq\> | Directory pattern for fastq files.                                   |
| --reference_genome    | \<hg19\>              | Path to the reference genome to which the samples will be mapped |
| --hla_ref | \<'hla_ref.fasta'>                      | Path to the hla reference genome |
| --hla_txt_file   | \<Allele_v2.txt\>                        | Path to the hla allele text file        |
