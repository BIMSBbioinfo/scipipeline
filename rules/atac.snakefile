from os.path import join
import os
import pysam

from utils.assemble_pseudogenome import create_pseudo_genome
from utils.split_reads import split_reads_by_barcode
from utils.split_reads import deduplicate_reads_by_barcode
from utils.count_matrix import sparse_count_reads_in_bins
from utils.count_matrix import sparse_count_reads_in_regions

# Trimmed reads
TRIM_PATTERN = join(OUT_DIR, 'atac_reads_trimmed')
FIRST_MATE_TRIMMED = TRIM_PATTERN + '_1.fastq'
SECOND_MATE_TRIMMED = TRIM_PATTERN + '_2.fastq'

# Reference genome and mapping index
MAPPING_RESULTS = join(OUT_DIR, '{reference}', 'atac.bam')
BAM_SORTED = join(OUT_DIR, '{reference}', 'atac_sorted.bam')

PSGENOME_OUTDIR = join(OUT_DIR, 'pseudogenome')



# ------------------------- #
# Adapter trimming
#
# Note: We set the -u parameter to 200
# to ensure that all reads are forwarded to the 
# next processing phase. Otherwise, matching
# the reads with the barcodes becomes a hassle
rule adapter_trimming:
    """Trim adapters using flexbar"""
    input: 
        first=config['samples']['FIRST_MATE'], 
        second=config['samples']['SECOND_MATE'], 
        adapters=config['adapters']
    output: 
        first=FIRST_MATE_TRIMMED, 
        second=SECOND_MATE_TRIMMED
    threads: 40
    log: join(LOG_DIR, 'flexbar.log')
    shell:
        "flexbar "
	+ "-r {input.first} -p {input.second} -t " + TRIM_PATTERN
        + " -a {input.adapters} "
        + " -f i1.8 -u 10 -ae RIGHT -at 1.0 --threads {threads} "
        + " --min-read-length 50 "
        + " > {log}"

# ------------------------- #
# Mapping to the reference genome

rule read_mapping:
    "Maps reads against reference genome"
    input: 
        first=FIRST_MATE_TRIMMED, 
        second=SECOND_MATE_TRIMMED
    output: MAPPING_RESULTS
    params: '/data/akalin/wkopp/bowtie2_indices/{reference}/{reference}'
    threads: 20
    log: join(LOG_DIR, '{reference}_bowtie2.log')
    shell:
        "bowtie2 -p {threads} -X 1500 --no-mixed " +
        "--no-discordant -x  {params} " +
        " -1 {input.first} -2 {input.second} 2> {log} | samtools view -bS - > " +
        "{output} "


# ------------------------- #
# Sort reads by name

rule sort_mapping_by_name:
    """Sort the reads by name"""
    input: MAPPING_RESULTS
    output: BAM_SORTED
    shell:
        "samtools sort -n {input} -o {output}"

# ------------------------- #
# Construct a pseudo genome in fasta format

rule make_pseudo_genomes:
    """Make pseudo genomes from indices"""
    input: 'meta/{barcode}.tab'
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.fasta')
    run:
        create_pseudo_genome(input[0], output[0])
       
# ------------------------- #
# Create bowtie2 index for pseudo genome
rule make_bowtie2_index_from_pseudo_genomes:
    """Make bowtie2 index from pseudo genomes"""
    params: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}')
    input: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.fasta')
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.1.bt2')
    shell:
        "bowtie2-build {input} {params}"
        

# ------------------------- #
# Mapping to pseudo genome

rule map_to_pseudo_genome:
    """Make bowtie2 index from pseudo genomes"""
    params: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}')
    input: 
        fastq = config['barcode_template'],
        index = join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.1.bt2')
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.bam')
    threads: 10
    log: join(LOG_DIR, 'pseudo_genome_{barcode}.log')
    shell:
        "bowtie2 -p {threads} -x {params} -U {input.fastq} 2> {log} | samtools view -bS - >" +
        " {output} "
        

# ------------------------- #
# Mapping to pseudo genome

rule sort_mapping_pseudo_genome_by_name:
    """Sort the reads by name"""
    input: join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}.bam')
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}_sorted.bam')
    shell:
        "samtools sort -n {input} -o {output}" 


# ------------------------- #
# Split the reads according to the barcodes
# 
rule split_reads_by_index:
    """Split reads by barcodes"""
    input:
       barcode_alns=expand(join(PSGENOME_OUTDIR, 
                                'pseudo_genome_{barcode}_sorted.bam'),
                                barcode=config['barcodes']),
       read_aln=BAM_SORTED
    output: join(OUT_DIR, "{reference}", "atac_barcoded.bam")
    params:
       min_mapq = 40,
       max_mismatches = 1
    run:
       split_reads_by_barcode(input.barcode_alns, 
                              input.read_aln, 
                              output[0],
                              params.min_mapq, params.max_mismatches)
       
#INPUT_ALL.append(expand(rules.split_reads_by_index.output, reference=config['reference']))
# ------------------------- #
# Sort the split reads

rule sort_split_reads:
    """Sort split reads"""
    input: join(OUT_DIR, "{reference}", "atac_barcoded.bam")
    output: join(OUT_DIR, "{reference}", "atac_barcoded_sorted.bam")
    shell: "samtools sort {input} -o {output}"

# ------------------------- #
# Deduplicate the split reads

rule deduplicate_split_reads:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "atac_barcoded_sorted.bam")
    output: join(OUT_DIR, "{reference}", "atac_barcode_dedup.bam")
    shell: "samtools rmdup {input} {output}"


# ------------------------- #
# Deduplicate the split reads by barcode

rule deduplicate_split_reads_by_barcode:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "atac_barcoded_sorted.bam")
    output: join(OUT_DIR, "{reference}", "atac_barcode_dedup_by_barcode.bam")
    run:
        deduplicate_reads_by_barcode(input[0], output[0])

INPUT_ALL.append(expand(rules.deduplicate_split_reads_by_barcode.output, reference=config['reference']))

# ------------------------- #
# Index the reads

rule index_deduplicate_split_reads:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "atac_barcode_dedup.bam")
    output: join(OUT_DIR, "{reference}", "atac_barcode_dedup.bam.bai")
    shell: "samtools index {input}"


# ------------------------- #
# report barcode frequencies

rule peak_calling_on_aggregate:
    input: join(OUT_DIR, "{reference}", "atac_barcode_dedup.bam")
    output: join(OUT_DIR, "{reference}", "macs2", "{reference}_peaks.narrowPeak"), join(OUT_DIR, "{reference}", "macs2", "{reference}_summits.bed")
    params: name='{reference}',
            outdir = join(OUT_DIR, "{reference}", "macs2")
    log: join(LOG_DIR, 'macs2_{reference}.log')
    shell: " macs2 callpeak --name {params.name} -t {input} -f BAMPE --nomodel --outdir {params.outdir} --call-summits --gsize dm 2> {log} "

INPUT_ALL.append(expand(rules.peak_calling_on_aggregate.output, reference=config['reference']))


# ------------------------- #
# Count reads in bins 

rule counting_reads_in_bins:
    """Counting reads per barcode"""
    input: 
        bams = join(OUT_DIR, "{reference}", "atac_barcode_dedup.bam"),
        bai = join(OUT_DIR, "{reference}", "atac_barcode_dedup.bam.bai")
    output: join(OUT_DIR, '{reference}', 'bin_counts_{binsize}.h5')
    log: join(LOG_DIR, 'count.bins{binsize}.{reference}.log')
    run: sparse_count_reads_in_bins(input.bams, int(wildcards.binsize), output[0])


INPUT_ALL.append(expand(rules.counting_reads_in_bins.output, reference=config['reference'], binsize=[2000, 5000, 10000]))

# ------------------------- #
# Count reads in regions

rule counting_reads_in_regions:
    """Counting reads per barcode"""
    input: 
        bams = join(OUT_DIR, "{reference}", "atac_barcode_dedup.bam"),
        bai = join(OUT_DIR, "{reference}", "atac_barcode_dedup.bam.bai"),
        regions = join(OUT_DIR, "{reference}", "macs2", "{reference}_summits.bed")
    output: join(OUT_DIR, '{reference}', 'region_counts.h5')
    params: flank=250
    log: join(LOG_DIR, 'count.region.{reference}.log')
    run: sparse_count_reads_in_regions(input.bams, input.regions, \
         output[0], flank=params.flank)


INPUT_ALL.append(expand(rules.counting_reads_in_regions.output, reference=config['reference']))

