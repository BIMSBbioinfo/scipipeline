# single-cell ATAC-seq

This repository contains code for analysing single-cell ATAC-seq dataset.

## Inputs to the pipeline

The input of the pipeline depends on whether combinatorial index was used
or whether the UMIs have been determined beforehand.

For combinatorial indexing the required input consists of

* Single- or paired-end reads for one or more samples in fastq format. The files may also be gzipped.
* Single- or paired-end reads of one or more barcoding rounds.
* One or more table or fasta files containing the reference barcodes.

Optional, adapter sequence in fasta format can be supplied for trimming.

In case barcodes have already been derived and associated with the reads
beforehand, e.g this might be the case if you have downloaded the dataset from an external source,
then the set of input files is slightly different

* Single- or paired-end reads for one or more samples in fastq format. The files may also be gzipped.
* (Optional) adapter sequences in fasta format.


## What the pipeline produces

The pipeline performs
* Read trimming (either using flexbar or trim_galore)
* Quality control (using fastqc before and after trimming)
* Read mapping (using bowtie2)
* Barcode correction using the reference barcodes (by aligning barcode reads against the reference barcodes with bowtie2)
* Barcode augmentation (appends the barcode to the RG tag)
* Barcode-aware deduplication (using Picard MarkDuplicates)
* Peak calling on the aggregated reads (using MACS2)
* Count matrix construction: Several count matrices are produced, 1) MACS2-peaks by cells and 2) genome-wide bins by cells. The binsizes and flanking windows can be adjusted in the worflowconfig.yaml file.

## Managing dependencies of the environment
We assume that the project root directory is located at `$PRJ_DIR`
which can be the root of the repository for example.

The necessary environment and tools to run the pipeline are defined
in guixr.manifest which is managed by guix.
To instantiate the environment the first time, run:

`guixr package --manifest=guixr.manifest --profile=$PRJ_DIR/.guix_profile`

Afterwards, the environment can be instantiated via
`source $PRJ_DIR/.guix_profile/etc/profile`.

Most tools are available within guix, except for Picard which is java-based.
To install Picard make sure that Java 1.8 is installed and proceed as follows:

```
git clone https://github.com/broadinstitute/picard.git
cd picard
./gradlew shadowJar
```
Afterwards edit `workflowconfig.yaml`
and add 
`picard_jarpath: <...>/picard/build/libs/picard.jar`


## Configure the pipeline

Currently the pipeline has three configuration points

* `workflowconfig.yaml`
* `samplesheet.tsv`
* `barcodesheet.tsv`

### workflowconfig.yaml
This file can be used to parametrize the workflow.
Below is a description of options:

```yaml
name: ProjectName

samples: "samplesheet.tsv"

# Specifies whether to perform sequencing error correction
# of the barcode reads or not (more details below).
# This part also allow to define minimum mapping quality
# and maximal number of mismatches for the barcode reads.
# The latter two options are irrelevant if correction: False
barcodes:
    sheet: "barcodesheet.tsv"
    correction: True/False
    min_mapq: 5
    max_mismatch: 2


# Specify whether the reads should be trimmmed with
# flexbar or with trim_galore
trim_reads: trim_galore/flexbar

# This options can be removed if no adapters are available
adapters: /path/to/adapters.fasta

# reference genomes to map against
# with the paths pointing to respective bowtie2 indices
# The pipeline can be used to align the reads to multiple reference
# genomes.
# The macs_gsize option defines the effective genome size for peak calling
# with macs2. removechroms allow to specify a list of chromosomes or
# patterns of chromosome names which should be removed for the downstream
# analysis.
min_mapq: 0
min_counts_per_barcode: 0
reference:
  hg38:
    bowtie2index: /data/akalin/wkopp/bowtie2_indices/hg38/genome
    macs_gsize: hs
    removechroms: [chrM, _random, chrUn]
  hg19:
    bowtie2index: /data/akalin/wkopp/bowtie2_indices/hg19/genome
    macs_gsize: hs
    removechroms: [chrM, _random, chrUn]


# binsizes to bin the reference genome
# in order to construct the count matrices
binsize:
  - 2000
  - 5000
  - 10000

# flank to extend the macs2 peak summits at both sides
# in order to construct the count matrix
peak_flank:
  - 250
  - 500

# output directory
output_dir: /path/outputdir

# this is the path to the java-based tool picard
picard_jarpath: /path/to/picard.jar
```
### samplesheet.tsv

This file contains a table 

|Name | read1 | read2 | barcodes |
|-----|-------|-------|----------|

* Name: Sample name
* read1: First mate of a paired-end sequencing run or a single-end sequenced file.
* read2: Second mate of a paired-end sequencing run or empty if it is a single-end sample.
* barcodes: ';' separated list of indices associated with the sample, e.g. `index1;index2;index3`. The index names are defined in barcodesheet.tsv.

### barcodesheet.tsv

The semantics of the barcodesheet.tsv differs when using barcode sequencing error correction or without correction.

In the case of sequencing error correction, which might be used for combinatorial indexing samples the file contains

| Name | read | reference |
|------|------|-----------|

* Name: Barcode or index name. For example, `index1`. These names must match with the index names listed in samplesheet.tsv.
* read: Fastq file containing the barcode reads.
* reference: Table containing one reference barcode per line or a fasta file containing the reference barcodes.

In case no barcode sequencing error correction should be performed, it is assumed that the reads have already
been associated with the barcodes. More specifically, 
the barcodes are assumed constitute a prefix of the read names in the samples.

The file barcodesheet.tsv then should contain a table 

| Readprefix | Name |
|------------|------|

* Readprefix: prefix string at the beginning of a read name which corresponds to the barcode
* Name: Alias name for the barcode.

## Run the snakemake pipeline

In order to run the pipeline issue the following command
```
snakemake -s main.snake
```

After the pipeline ran through additionally run

```
snakemake -s main.snake --report results.html
```

The latter command will create a html report containing some result figures of the analysis,
including barcode frequencies, fragment length distribution, etc.

## Results

The results are organized in the specified output directory as follows:

```
|-- trimmed_1.fq, trimmed_2.fq
|-- <genome_name>
|     |-- sample.bam
|     |-- sample.barcoded.minmapq{XX}.dedup.mincount{XX}.bam
|     |-- sample.barcoded.minmapq{XX}.dedup.mincount{XX}.bw
|     |-- countmatrix
|     |      |-- *.tab: count matrix as sparse dataset
|     |      |-- *.bed: regions associated with the count matrix
|     |      \-- *.bed: regions associated with the count matrix
|     |
|     |-- macs2: macs2 results
|     \-- report: summary figures and tables
|
|-- logs
|-- fastqc_raw
|-- fastqc_trimmed
|-- multiqc_data
|-- multiqc_report.html
```

## Hints

### Whether to use flexbar or trim_galore
 We suggest to use flexbar with a set of known adapter sequences.
 If the adapters are not known, we suggest to use trim_galore. The latter tool
 is able to infer commonly used sequencing adapters.

