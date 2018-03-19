from os.path import join

from utils.assemble_pseudogenome import create_pseudo_genome
from utils.split_reads import split_reads_by_barcode
from utils.split_reads import obtain_barcode_frequencies

# Trimmed reads
TRIM_PATTERN = join(OUT_DIR, 'atac_reads_trimmed')
FIRST_MATE_TRIMMED = TRIM_PATTERN + '_1.fastq'
SECOND_MATE_TRIMMED = TRIM_PATTERN + '_2.fastq'

# Reference genome and mapping index
FLY_BAM = join(OUT_DIR, 'atac.bam')
FLY_BAM_SORTED = join(OUT_DIR, 'atac_sorted.bam')
FLY_BOWTIE2_INDEX = '/data/ohler/Mahmoud/info-dm6/bowtie2/dm6'

PSGENOME_OUTDIR = join(OUT_DIR, 'pseudogenome')
SPLIT_OUTPUT_DIR = join(OUT_DIR, 'splitreads')



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
        first=FIRST_MATE, 
        second=SECOND_MATE,
        adapters=ADAPTERS
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

INPUT_ALL.append(rules.adapter_trimming.output)

# ------------------------- #
# Mapping to fly genome

rule read_mapping_fly:
    input: 
        first=FIRST_MATE_TRIMMED, 
        second=SECOND_MATE_TRIMMED
    output: FLY_BAM
    threads: 20
    log: join(LOG_DIR, 'fly_bowtie2.log')
    shell:
        "bowtie2 -p {threads} -X 1500 --no-mixed " +
        "--no-discordant -x " + FLY_BOWTIE2_INDEX +
        " -1 {input.first} -2 {input.second} 2> {log} | samtools view -bS - > " +
        "{output} "


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

rule make_pseudo_genomes:
    """Make pseudo genomes from indices"""
    output: [join(PSGENOME_OUTDIR, 'pseudo_genome_I{}.fasta'.format(sample)) for sample in [1,2,3,4]]
    run:
        create_pseudo_genome(PSGENOME_OUTDIR)
       
INPUT_ALL.append(rules.make_pseudo_genomes.output)

# ------------------------- #
# Create bowtie2 index for pseudo genome
rule make_bowtie2_index_from_pseudo_genomes:
    """Make bowtie2 index from pseudo genomes"""
    params: join(PSGENOME_OUTDIR, 'pseudo_genome_I{sample}')
    input: join(PSGENOME_OUTDIR, 'pseudo_genome_I{sample}.fasta')
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_I{sample}.1.bt2')
    shell:
        "bowtie2-build {input} {params}"
        
INPUT_ALL.append(expand(join(PSGENOME_OUTDIR, 'pseudo_genome_I{sample}.1.bt2'), sample=[1,2,3,4]))


# ------------------------- #
# Mapping to pseudo genome

rule map_to_pseudo_genomes:
    """Make bowtie2 index from pseudo genomes"""
    params: join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}')
    input: 
        fastq = BARCODES
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}.bam')
    threads: 10
    log: join(LOG_DIR, 'pseudo_genome_{sample}.log')
    shell:
        "bowtie2 -p {threads} -x {params} -U {input.fastq} 2> {log} | samtools view -bS - >" +
        " {output} "
        
INPUT_ALL.append(expand(join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}.bam'), 
                        sample=['I1', 'I2', 'I3', 'I4']))

# ------------------------- #
# Mapping to pseudo genome

rule sort_mapping_pseudogenome:
    """Sort the reads by name"""
    input: join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}.bam')
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}_sorted.bam')
    shell:
        "samtools sort -n {input} -o {output}" 

INPUT_ALL.append(expand(join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}_sorted.bam'), 
                        sample=['I1', 'I2', 'I3', 'I4']))

# ------------------------- #
# Sort the split reads

rule sort_split_reads:
  input: join(SPLIT_OUTPUT_DIR, "{barcode}.bam")
  output: join(SPLIT_OUTPUT_DIR + '_sorted', "{barcode}.bam")
  shell: "samtools sort {input} -o {output}"

# ------------------------- #
# Deduplicate the split reads

rule deduplicate_split_reads:
  input: join(SPLIT_OUTPUT_DIR + '_sorted', "{barcode}.bam")
  output: join(SPLIT_OUTPUT_DIR + '_deduplicated', "{barcode}.bam")
  shell: "samtools rmdup {input} {output}"

#INPUT_ALL.append(dynamic(join(SPLIT_OUTPUT_DIR + '_deduplicated', '{barcode}.bam')))
#INPUT_ALL.append(dynamic(join(SPLIT_OUTPUT_DIR, '{barcode}.bam')))

# ------------------------- #
# report barcode frequencies

rule report_barcode_frequencies:
    input: 
        original = dynamic(join(SPLIT_OUTPUT_DIR, "{barcode}.bam")),
        deduplicated = dynamic(join(SPLIT_OUTPUT_DIR + '_deduplicated', "{barcode}.bam"))
    output:
        join(OUT_DIR,"barcode_frequencies.tab")
    run:
        obtain_barcode_frequencies(input.original, input.deduplicated, output[0])
        
INPUT_ALL.append(rules.report_barcode_frequencies.output)
# ------------------------- #
# Split the reads according to the barcodes
# 
rule split_reads_by_index:
    input:
       barcode_alns=expand(join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}_sorted.bam'),
                        sample=['I1', 'I2', 'I3', 'I4']),
       read_aln=FLY_BAM_SORTED
    output: dynamic(join(SPLIT_OUTPUT_DIR, "{barcode}.bam"))
    run:
       split_reads_by_barcode(input.barcode_alns, 
                              input.read_aln, 
                              os.path.dirname(output[0]), 1000)
       
INPUT_ALL.append(rules.split_reads_by_index.output)

