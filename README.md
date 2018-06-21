# single-cell ATAC-seq

This repository contains code for analysing single-cell ATAC-seq dataset
based on a combinatorial indexing approach.

## What the pipeline produces

The pipeline takes as input 

1. single- or paired-end reads for one or more samples in fastq format.
1a. one or more reference genomes to map reads against
1b. adapter sequences in fasta format to trim the reads
2. one or more sets of barcode/UMI reads in fastq format
2b. reference barcodes to correct the sequenced barcodes as table or fasta file.

The pipeline performs
* reads mapping
* barcode correction using the reference barcodes
* read splitting by barcode
* deduplication of alignments by barcode
* peak calling on the aggregated reads across barcodes
* computation of count matrices: peaks by cells and bins by cells bins by cells allows to bin the genome in equally long windows.


## Setup the environment
We assume that the project root directory is located at `$PRJ_DIR`.

The environment is managed by guix.
To instantiate the environment the first time, run:

`guixr package --manifest=guixr.manifest --profile=$PRJ_DIR/.guix_profile`

Afterwards, the environment can be instantiated via
`source $PRJ_DIR/.guix_profile/etc/profile`.

## Run the snakemake pipeline

The dataset is processed using a snakemake pipeline. To reproduce the results
locally run
```
snakemake -s main.snake
```

Alternatively, the pipeline can be launched on the cluster using
```
snakemake -s main.snake -j <njobs> -V -cwd
```

To run the analysis, `ssh max-cluster2` and
afterwards `qrsh`

The pipeline rules are specified in the `rules` folder,
`src` contains utilities and functions that are used in the
analysis or in the pipeline.
`analysis` contains scripts that result in reports or plots.
For example, to address project specific questions
like hypothesis tests or exploratory analysis or other

