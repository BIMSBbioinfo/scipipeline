
# ------------------------- #
# Sort reads by name

rule sort_mapping_by_name:
    """Sort the reads by name"""
    input: join(OUT_DIR, '{reference}', '{sample}.cleanchrom.bam')
    output: join(OUT_DIR, '{reference}', '{sample}.namesorted.bam')
    resources:
        mem_mb=10000
    wildcard_constraints:
        sample="[\w]+"
    shell:
        "samtools sort -n {input} -o {output}"

# ------------------------- #
# Construct a pseudo genome in fasta format

rule make_pseudo_genomes:
    """Make pseudo genomes from indices"""
    input: lambda wildcards: barcodes[barcodes.Name == wildcards.barcode].reference.tolist()
    output: join(PSGENOME_OUTDIR, '{barcode}.fasta')
    resources:
      mem_mb=100
    wildcard_constraints:
        barcode="[\w\d]+"
    run:
        create_pseudo_genome(input[0], output[0])

# ------------------------- #
# Create bowtie2 index for pseudo genome
rule make_bowtie2_index_from_pseudo_genomes:
    """Make bowtie2 index from pseudo genomes"""
    params: join(PSGENOME_OUTDIR, '{barcode}')
    input: join(PSGENOME_OUTDIR, '{barcode}.fasta')
    output: join(PSGENOME_OUTDIR, '{barcode}.1.bt2')
    resources:
        mem_mb=100
    shell:
        "bowtie2-build {input} {params}"

# ------------------------- #
# Mapping to pseudo genome
def _bowtie_input_type_bc(wildcards):
    filename = barcodes[barcodes.Name==wildcards.barcode].read.tolist()[0]
    return bowtie_input_filetype_option(filename)

rule map_to_pseudo_genome:
    """Make bowtie2 index from pseudo genomes"""
    params:
        genome = join(PSGENOME_OUTDIR, '{barcode}'),
        filetype = _bowtie_input_type_bc
    input:
        fastq = lambda wildcards: barcodes[barcodes.Name==wildcards.barcode].read.tolist(),
        index = join(PSGENOME_OUTDIR, '{barcode}.1.bt2')
    output: join(PSGENOME_OUTDIR, 'barcode.{barcode}.bam')
    threads: 10
    resources:
       mem_mb=500
    log: join(LOG_DIR, '{barcode}.log')
    wildcard_constraints:
        barcode="[\w\d]+"
    shell:
        "bowtie2 -p {threads} -x {params.genome} " +
        "{params.filetype} -U {input.fastq} 2> {log} " +
        " | samtools view -bS - > {output} "


# ------------------------- #
# Mapping to pseudo genome

rule sort_mapping_pseudo_genome_by_name:
    """Sort the reads by name"""
    input: join(PSGENOME_OUTDIR, 'barcode.{barcode}.bam')
    output: join(PSGENOME_OUTDIR, 'barcode.{barcode}.namesorted.bam')
    resources:
        mem_mb=10000
    threads: 2
    wildcard_constraints:
        barcode="[\w\d]+"
    shell:
        "samtools sort -n {input} -o {output}"


# ------------------------- #
# Split the reads according to the barcodes
#

rule split_reads_by_index:
    """Split reads by barcodes"""
    input:
       barcode_alns= lambda wc:
                        expand(join(PSGENOME_OUTDIR,
                                'barcode.{barcode}.namesorted.bam'),
                                barcode=samples[samples.Name==wc.sample].barcodes.tolist()[0].split(';')),
       read_aln=join(OUT_DIR, '{reference}', '{sample}.namesorted.bam')
    output: join(OUT_DIR, "{reference}", "{sample}.barcoded.bam")
    resources:
        mem_mb=1000
    threads: 2
    params:
       min_mapq = config['barcodes']['min_mapq'],
       max_mismatches = config['barcodes']['max_mismatch']
    run:
       split_reads_by_barcode(input.barcode_alns,
                              input.read_aln,
                              output[0],
                              params.min_mapq, params.max_mismatches)
