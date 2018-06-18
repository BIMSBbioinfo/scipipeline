from os.path import join
import os
import pysam

from utils.assemble_pseudogenome import create_pseudo_genome
from utils.split_reads import split_reads_by_barcode
from utils.split_reads import deduplicate_reads
from utils.count_matrix import count_reads_in_bins
from utils.count_matrix import count_reads_in_regions

# Reference genome and mapping index
MAPPING_RESULTS = join(OUT_DIR, '{reference}', '{sample}.bam')
BAM_SORTED = join(OUT_DIR, '{reference}', '{sample}.sorted.bam')

PSGENOME_OUTDIR = join(OUT_DIR, 'barcodes')


# ------------------------- #
# Adapter trimming
#
# Note: We set the -u parameter to 200
# to ensure that all reads are forwarded to the 
# next processing phase. Otherwise, matching
# the reads with the barcodes becomes a hassle

# Trimmed reads
TRIM_PATTERN = join(OUT_DIR, '{wildcards.sample}_trimmed')

if config['pairedend']:
    FIRST_MATE_TRIMMED = join(OUT_DIR, '{sample}_trimmed_1.fastq')
    SECOND_MATE_TRIMMED = join(OUT_DIR, '{sample}_trimmed_2.fastq')

    rule adapter_trimming:
        """Trim adapters using flexbar"""
        input: 
            first=lambda wildcards: config['samples'][wildcards.sample]['read1'], 
            second=lambda wildcards: config['samples'][wildcards.sample]['read2'], 
            adapters=config['adapters']
        output: [FIRST_MATE_TRIMMED, SECOND_MATE_TRIMMED]
        threads: 40
        log: join(LOG_DIR, 'flexbar_{sample}.log')
        shell:
            "flexbar "
	    + "-r {input.first} -p {input.second} -t "+TRIM_PATTERN
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
        log: join(LOG_DIR, '{sample}_{reference}_bowtie2.log')
        wildcard_constraints:
            sample="\w+"
        shell:
            "bowtie2 -p {threads} -X 1500 --no-mixed " +
            "--no-discordant -x  {params} " +
            " -1 {input.first} -2 {input.second} 2> {log} |" +
            " samtools view -bS - > {output} "

else:
    FIRST_MATE_TRIMMED = join(OUT_DIR, '{sample}_trimmed.fastq')

    rule adapter_trimming:
        """Trim adapters using flexbar"""
        input: 
            first=lambda wildcards: config['samples'][wildcards.sample]['read2'], 
            adapters=config['adapters']
        output: [FIRST_MATE_TRIMMED]
        threads: 40
        log: join(LOG_DIR, 'flexbar.log')
        shell:
            "flexbar "
	    + "-r {input.first} -t " + TRIM_PATTERN
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
        output: MAPPING_RESULTS
        params: '/data/akalin/wkopp/bowtie2_indices/{reference}/{reference}'
        threads: 20
        log: join(LOG_DIR, '{reference}_bowtie2.log')
        wildcard_constraints:
            sample="\w+"
        shell:
            "bowtie2 -p {threads} -X 1500 --no-mixed " +
            "--no-discordant -x  {params} " +
            " -U {input.first} 2> {log} | " +
            "samtools view -bS - > {output} "


# ------------------------- #
# Sort reads by name


rule sort_mapping_by_name:
    """Sort the reads by name"""
    input: MAPPING_RESULTS
    output: BAM_SORTED
    wildcard_constraints:
        sample="[\w]+"
    shell:
        "samtools sort -n {input} -o {output}"

# ------------------------- #
# Construct a pseudo genome in fasta format

rule make_pseudo_genomes:
    """Make pseudo genomes from indices"""
    input: lambda wildcards: config['barcodes'][wildcards.barcode]['reference']
    output: join(PSGENOME_OUTDIR, '{barcode}.fasta')
    run:
        create_pseudo_genome(input[0], output[0])
       
# ------------------------- #
# Create bowtie2 index for pseudo genome
rule make_bowtie2_index_from_pseudo_genomes:
    """Make bowtie2 index from pseudo genomes"""
    params: join(PSGENOME_OUTDIR, '{barcode}')
    input: join(PSGENOME_OUTDIR, '{barcode}.fasta')
    output: join(PSGENOME_OUTDIR, '{barcode}.1.bt2')
    shell:
        "bowtie2-build {input} {params}"
        

# ------------------------- #
# Mapping to pseudo genome

rule map_to_pseudo_genome:
    """Make bowtie2 index from pseudo genomes"""
    params: join(PSGENOME_OUTDIR, '{barcode}')
    input: 
        fastq = lambda wildcards: config['barcodes'][wildcards.barcode]['reads'],
        index = join(PSGENOME_OUTDIR, '{barcode}.1.bt2')
    output: join(PSGENOME_OUTDIR, 'barcode.{barcode}.bam')
    threads: 10
    log: join(LOG_DIR, '{barcode}.log')
    shell:
        "bowtie2 -p {threads} -x {params} -U {input.fastq} 2> {log} " +
        " | samtools view -bS - > {output} "
        

# ------------------------- #
# Mapping to pseudo genome

rule sort_mapping_pseudo_genome_by_name:
    """Sort the reads by name"""
    input: join(PSGENOME_OUTDIR, 'barcode.{barcode}.bam')
    output: join(PSGENOME_OUTDIR, 'barcode.{barcode}.sorted.bam')
    wildcard_constraints:
        barcode="[\w\d]+"
    shell:
        "samtools sort -n {input} -o {output}" 


# ------------------------- #
# Split the reads according to the barcodes
# 
print(expand(join(PSGENOME_OUTDIR,
                                'barcode.{barcode}.sorted.bam'),
                                barcode=config['barcodes'].keys()))
print(BAM_SORTED)
rule split_reads_by_index:
    """Split reads by barcodes"""
    input:
       barcode_alns=expand(join(PSGENOME_OUTDIR, 
                                'barcode.{barcode}.sorted.bam'),
                                barcode=config['barcodes'].keys()),
       read_aln=BAM_SORTED
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.bam")
    params:
       min_mapq = 40,
       max_mismatches = 1
    run:
       split_reads_by_barcode(input.barcode_alns, 
                              input.read_aln, 
                              output[0],
                              params.min_mapq, params.max_mismatches)
       
# ------------------------- #
# Sort the split reads

rule sort_split_reads:
    """Sort split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.sorted.bam")
    shell: "samtools sort {input} -o {output}"

# ------------------------- #
# Deduplicate the split reads by barcode

rule deduplicate_split_reads_by_barcode:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.sorted.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam")
    run:
        deduplicate_reads(input[0], output[0])

INPUT_ALL.append(expand(rules.deduplicate_split_reads_by_barcode.output, reference=config['reference'], sample=config['samples'].keys()))

# ------------------------- #
# Index the reads

rule index_deduplicate_split_reads:
    """Deduplicate split reads"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam")
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam.bai")
    shell: "samtools index {input}"


# ------------------------- #
# report barcode frequencies

rule peak_calling_on_aggregate:
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam")
    output: join(OUT_DIR, "{reference}", "macs2", "{sample}_peaks.narrowPeak"), join(OUT_DIR, "{reference}", "macs2", "{sample}_summits.bed")
    params: name='{reference}',
            outdir = join(OUT_DIR, "{reference}", "macs2")
    log: join(LOG_DIR, 'macs2_{sample}.log')
    shell: 
      " macs2 callpeak --name {params.name} -t {input} -f " +
      "BAM" + "PE" if config['pairedend'] else "" +
      " --nomodel --outdir {params.outdir} --call-summits --gsize dm 2> {log} "

INPUT_ALL.append(expand(rules.peak_calling_on_aggregate.output, reference=config['reference'], sample=config['samples'].keys()))


# ------------------------- #
# Count reads in bins 

rule counting_reads_in_bins:
    """Counting reads per barcode"""
    input: 
        bams = join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam"),
        bai = join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam.bai")
    output: join(OUT_DIR, '{reference}', 'countmatrix', '{sample}_bin_counts_{binsize}.h5')
    run: count_reads_in_bins(input.bams, int(wildcards.binsize), output[0])


INPUT_ALL.append(expand(rules.counting_reads_in_bins.output, 
                        reference=config['reference'], 
                        binsize=[2000, 5000, 10000], 
                        sample=config['samples'].keys()))

# ------------------------- #
# Count reads in regions

rule counting_reads_in_regions:
    """Counting reads per barcode"""
    input: 
        bams = join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam"),
        bai = join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam.bai"),
        regions = join(OUT_DIR, "{reference}", "macs2", "{sample}_summits.bed")
    output: join(OUT_DIR, '{reference}', 'countmatrix', '{sample}_peak_counts.h5')
    params: flank=250
    run: count_reads_in_regions(input.bams, input.regions, \
         output[0], flank=params.flank)


INPUT_ALL.append(expand(rules.counting_reads_in_regions.output, reference=config['reference'], sample=config['samples'].keys()))

