# single-cell combinatorial indexing ATAC-seq pipeline

<div align="center">
<img src="drawing.png" alt="sciATACseq Pipeline" width=70%></>
</div>


## Hallmarks of the pipeline
The pipeline is implemented using snakemake and has the following features:

* Compatible with single-/ paired-end sequencing data
* Trimming with / without known adapters
* Alignment against one or more genomes
* Barcode correction using reference barcodes
* Supports One or more combinatorial indexing rounds
* Supports multicore or grid-engine environment


## Input description

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


## What the pipeline does

The pipeline performs
* Read trimming (either using flexbar or trim_galore depending whether adapter sequences have been provided)
* Quality control (using fastqc before and after trimming)
* Read mapping (using bowtie2)
* Barcode correction using the reference barcodes (by aligning barcode reads against the reference barcodes with bowtie2)
* Barcode augmentation (appends the barcode ID to the RG tag of the mapped reads)
* Barcode-aware deduplication (using Picard MarkDuplicates)
* Peak calling on the aggregated reads (using MACS2)
* Count matrix construction: Several count matrices are produced, 1) MACS2-peaks by cells and 2) genome-wide bins by cells. The binsizes and flanking windows can be adjusted in the worflowconfig.yaml file.
* Determines barcode statistics
* Generates MultiQC report (in HTML format)
* Generate snakemake report (in HTML format)

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

The pipeline has three configuration points

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


# If the adapters are known, put them in a fasta file under the adapters field.
# flexbar will make use of them.
# On the other hand, if the adapters are not know, remove the adapters field.
# In this case, trim_galore will try to infer the adapter sequences.
adapters: /data/ohler/Scott/Asli_Scripts/sciAtacAdapters.fa

# The next block specifies the
# reference genomes to map against
# with the paths pointing to respective bowtie2 indices.
# The pipeline can be used to align the reads to multiple reference
# genomes.
# The macs_gsize option defines the effective genome size for peak calling
# with macs2. removechroms allow to specify a list of chromosomes or
# patterns of chromosome names which should be removed for the downstream
# analysis.
# Finally, under the annotation section, extra genomic regions can be specified
# for which countmatrices should be generated.
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
    annotation:
      markergenes: path/to/markergenes.bed
      markerenhancers: path/to/markerenhancers.bed

# Indictes that both ends of a paired-end alignment should be counted
# separatly for the construction of the count matrix.
# If True, both mates are counted at the 5' end. Otherwise,
# the midpoint between the two mates is counted once.
count_both_ends: True

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

The pipeline can be run locally in a grid-engine cluster environment.

In order to run the pipeline locally, use
```
snakemake -s main.snake --cores <ncores>
```

Alternatively, the pipeline can be run on a cluster environment.
To this end, the pipeline contains launcher scripts.

If you are allowed to run the pipeline directly from the head node
invoke `bash runpipeline_local.bash`. In that case, snakemake will run on
the current node and it will submit the individual jobs that need to be carried
out via `qsub` to the cluster.

If you are not allowed to run the pipeline on the head node
you can invoke `bash runpipeline_remote.bash`. This will start snakemake itself via qsub
on the cluster. From there `runpipeline_local.bash` will be used.
Note that this assumes that `qsub` can be issues from the node snakemake has been started.

After the pipeline has completed a report document can be generated via

```
snakemake -s main.snake --report results.html
```

The latter command will create a html report containing summary figures of the analysis,
including barcode frequencies, fragment length distribution, etc.

## Results

The results are organized in the specified output directory as follows:

```
|-- logs: outputs, warnings, errors of individual tools
|-- multiqc_data
|-- multiqc_report.html
|-- barcodes: barcode alignments
|-- <sample>
|    |-- fastqc_<raw/trimmed>: quality control
|    |-- trimmed: trimmed reads
|    |-- <reference>
|    |    |-- mapping: alignment files
|    |    |-- peaks: peak calling results
|    |    |-- countmatrix
|    |    |   |-- *.tab: count matrix as sparse dataset
|    |    |   |-- *.bed: regions associated with the count matrix
|    |    |   \-- *.bed: regions associated with the count matrix
|    |    |-- report: summary figures and tables
```

