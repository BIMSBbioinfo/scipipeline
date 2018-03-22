# single-cell ATAC-seq

This repository contains code for analysing single-cell ATAC-seq dataset
based on a combinatorial indexing approach.

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

