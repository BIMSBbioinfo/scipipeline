# single-cell ATAC-seq

This repository contains code for analysing single-cell ATAC-seq dataset
based on a combinatorial indexing approach.

## Setup the environment
The environment is managed by guix.
To instantiate the environment run:

`guixr package --manifest=guixr.manifest --profile=./.guix_profile`

## Run the snakemake pipeline

The dataset is processed using a snakemake pipeline. To reproduce the results
invoke:
```
snakemake -s main.snake
```

To run the analysis, `ssh max-cluster2` and
afterwards `qrsh`

The pipeline rules are specified in the `rules` folder,
`src` contains utilities and functions that are used in the
analysis or in the pipeline.
`analysis` contains scripts that result in reports or plots.
For example, to address project specific questions
like hypothesis tests or exploratory analysis or other

##
