import os

# Sequencing adapters for trimming
ADAPTERS = '/data/ohler/Scott/Asli_Scripts/sciAtacAdapters.fa'

# Definition of input and output files:
# ATAC-seq reads
FIRST_MATE = '/scratch/AG_Ohler/Scott/sciATACtest/noSampleSheet/Undetermined_S0_R1_001.fastq.gz'
SECOND_MATE = '/scratch/AG_Ohler/Scott/sciATACtest/noSampleSheet/Undetermined_S0_R2_001.fastq.gz'

# Trimmed reads
TRIM_PATTERN = OUT_DIR + 'atac_reads_trimmed'
FIRST_MATE_TRIMMED = TRIM_PATTERN + '_1.fastq'
SECOND_MATE_TRIMMED = TRIM_PATTERN + '_2.fastq'

# Reference genome and mapping index
FLY_BAM = OUT_DIR + 'atac.bam'
FLY_BAM_SORTED = OUT_DIR + 'atac_sorted.bam'
FLY_BOWTIE2_INDEX = '/data/ohler/Mahmoud/info-dm6/bowtie2/dm6'

PSEUDOGENOME_MAPPING = ""
PSEUDOGENOME = ""

# ------------------------- #
# Adapter trimming

rule adapter_trimming:
    """Trim adapters using flexbar"""
    input: 
        first=FIRST_MATE, 
        second=SECOND_MATE
    output: 
        first=FIRST_MATE_TRIMMED, 
        second=SECOND_MATE_TRIMMED
    shell:
        "flexbar "
	+ "-r {input.first} -p {input.second} -t " + TRIM_PATTERN
        + " -f i1.8 -u 10 -ae RIGHT -at 1.0"

INPUT_ALL.append(rules.adapter_trimming.output)

# ------------------------- #
# Mapping to fly genome

rule read_mapping_fly:
    input: 
        first=FIRST_MATE_TRIMMED, 
        second=SECOND_MATE_TRIMMED
    output: FLY_BAM
    threads: 5
    shell:
        "bowtie2 -p {threads} -X 1500 --no-mixed " +
        "--no-discordant -x " + FLY_BOWTIE2_INDEX +
        " -1 {input.first} -2 {input.second} | samtools view -bS - > " +
        "{output}"


INPUT_ALL.append(rules.read_mapping_fly.output)

# ------------------------- #
# Mapping to fly genome

rule sort_mapping_fly:
    """Sort the reads by name"""
    input: FLY_BAM
    output: FLY_BAM_SORTED
    shell:
        "samtools sort -n {input} -o {output}"

INPUT_ALL.append(rules.sort_mapping_fly.output)

# ------------------------- #
# Construct a pseudo genome in fasta format

# ------------------------- #
# Create bowtie2 index for pseudo genome

# ------------------------- #
# Mapping to pseudo genome


