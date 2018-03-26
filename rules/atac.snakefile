from os.path import join

from utils.assemble_pseudogenome import create_pseudo_genome
from utils.split_reads import split_reads_by_barcode
from utils.split_reads import obtain_barcode_frequencies
from utils.split_reads import plot_barcode_frequencies
from utils.split_reads import species_specificity

# Trimmed reads
TRIM_PATTERN = join(OUT_DIR, 'atac_reads_trimmed')
FIRST_MATE_TRIMMED = TRIM_PATTERN + '_1.fastq'
SECOND_MATE_TRIMMED = TRIM_PATTERN + '_2.fastq'

# Reference genome and mapping index
MAPPING_RESULTS = join(OUT_DIR, '{reference}', 'atac.bam')
BAM_SORTED = join(OUT_DIR, '{reference}', 'atac_sorted.bam')

PSGENOME_OUTDIR = join(OUT_DIR, 'pseudogenome')
SPLIT_OUTPUT_DIR = join(OUT_DIR, '{reference}', 'splitreads')



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

#INPUT_ALL.append(rules.adapter_trimming.output)

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


#INPUT_ALL.append(expand(rules.read_mapping.output, reference=config['reference']))

# ------------------------- #
# Sort reads by name

rule sort_mapping_by_name:
    """Sort the reads by name"""
    input: MAPPING_RESULTS
    output: BAM_SORTED
    shell:
        "samtools sort -n {input} -o {output}"

#INPUT_ALL.append(expand(rules.sort_mapping_by_name.output, reference=config['reference']))

# ------------------------- #
# Construct a pseudo genome in fasta format

rule make_pseudo_genomes:
    """Make pseudo genomes from indices"""
    input: 'meta/{barcode}.tab'
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.fasta')
    run:
        create_pseudo_genome(input[0], output[0])
       
#INPUT_ALL.append(expand(rules.make_pseudo_genomes.output, barcode=config['barcodes']))

# ------------------------- #
# Create bowtie2 index for pseudo genome
rule make_bowtie2_index_from_pseudo_genomes:
    """Make bowtie2 index from pseudo genomes"""
    params: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}')
    input: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.fasta')
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.1.bt2')
    shell:
        "bowtie2-build {input} {params}"
        
#INPUT_ALL.append(expand(join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.1.bt2'), barcode=config['barcodes']))


# ------------------------- #
# Mapping to pseudo genome

rule map_to_pseudo_genome:
    """Make bowtie2 index from pseudo genomes"""
    params: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}')
    input: 
        fastq = BARCODE_TEMPLATE,
        index = join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.1.bt2')
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.bam')
    threads: 10
    log: join(LOG_DIR, 'pseudo_genome_{barcode}.log')
    shell:
        "bowtie2 -p {threads} -x {params} -U {input.fastq} 2> {log} | samtools view -bS - >" +
        " {output} "
        
#INPUT_ALL.append(expand(join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}.bam'), 
                        #barcode=['I1', 'I2', 'I3', 'I4']))

# ------------------------- #
# Mapping to pseudo genome

rule sort_mapping_pseudo_genome_by_name:
    """Sort the reads by name"""
    input: join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}.bam')
    output: join(PSGENOME_OUTDIR, 'pseudo_genome_{sample}_sorted.bam')
    shell:
        "samtools sort -n {input} -o {output}" 

#INPUT_ALL.append(expand(join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}_sorted.bam'), 
                        #barcode=config['barcodes']))

# ------------------------- #
# Sort the split reads

rule sort_split_reads:
    """Sort split reads"""
    input: join(SPLIT_OUTPUT_DIR, "{combar}.bam")
    output: join(SPLIT_OUTPUT_DIR + '_sorted', "{combar}.bam")
    shell: "samtools sort {input} -o {output}"

# ------------------------- #
# Deduplicate the split reads

rule deduplicate_split_reads:
    """Deduplicate split reads"""
    input: join(SPLIT_OUTPUT_DIR + '_sorted', "{combar}.bam")
    output: join(SPLIT_OUTPUT_DIR + '_deduplicated', "{combar}.bam")
    shell: "samtools rmdup {input} {output}"


# ------------------------- #
# report barcode frequencies

rule report_barcode_frequencies:
    """Report barcode frequencies"""
    input: 
        original = dynamic(join(SPLIT_OUTPUT_DIR, "{combar}.bam")),
        deduplicated = dynamic(join(SPLIT_OUTPUT_DIR + '_deduplicated', "{combar}.bam"))
    output:
        join(OUT_DIR, "report", "barcode_frequencies.{reference}.tab")
    run:
        obtain_barcode_frequencies(input.original, input.deduplicated, output[0])
        
INPUT_ALL.append(expand(rules.report_barcode_frequencies.output, reference=config['reference']))

rule plot_barcode_freqs:
    """Plot barcode frequency"""
    input: join(OUT_DIR, "report", "barcode_frequencies.{reference}.tab")
    output: join(OUT_DIR, "report", "barcode_frequencies_{reference}.png")
    run:
        plot_barcode_frequencies(input[0], output[0])

INPUT_ALL.append(expand(rules.plot_barcode_freqs.output, reference=config['reference']))

# ------------------------- #
# Split the reads according to the barcodes
# 
rule split_reads_by_index:
    """Split reads by barcodes"""
    input:
       barcode_alns=expand(join(PSGENOME_OUTDIR, 'pseudo_genome_{barcode}_sorted.bam'),
                        barcode=config['barcodes']),
       read_aln=BAM_SORTED
    output: dynamic(join(SPLIT_OUTPUT_DIR, "{combar}.bam"))
    params:
       max_open_files = 1000,
       min_mapq = 40,
       max_mismatches = 1
    run:
       split_reads_by_barcode(input.barcode_alns, 
                              input.read_aln, 
                              os.path.dirname(output[0]), params.max_open_files,
                              params.min_mapq, params.max_mismatches)
       
#INPUT_ALL.append(dynamic(expand(join(SPLIT_OUTPUT_DIR, "{combar}.bam"), reference=config['reference'])))

rule read_cooccurrence_by_species:
    """How frequently do reads map uniquely to one species or to both?"""
    input: expand(BAM_SORTED, reference=config['reference'])
    output: expand(join(OUT_DIR, "report", "read_cooccurrence.{suffix}"), suffix=['tab', 'png'])
    run:
        species_specificity(input[0], input[1], output[0].split('.')[0], expand('{reference}', reference=config['reference']))

INPUT_ALL.append(rules.read_cooccurrence_by_species.output)
