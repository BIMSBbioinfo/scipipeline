from utils.cleanup_alignments import remove_chroms
from utils.cleanup_alignments import make_genome_size

def _bowtie_input_type_read(wildcards):
    filename = samples[samples.Name == wildcards.sample].read1.tolist()[0]
    return bowtie_input_filetype_option(filename)

rule read_mapping:
    "Maps reads against reference genome"
    input: get_mapping_inputs
    output: join(OUT_DIR, '{sample}', '{reference}', 'mapping', 'sample.bam')
    wildcard_constraints: sample='\w+'
    params:
        genome=lambda wildcards: config['reference'][wildcards.reference]['bowtie2index'],
        paired=is_paired,
        filetype=_bowtie_input_type_read
    threads: 20
    resources:
       mem_mb=1000
    log: join(LOG_DIR, '{sample}_{reference}_bowtie2.log')
    run:
        cmd = 'bowtie2'
        cmd = "bowtie2 -p {threads} -X 2000 --no-mixed --no-discordant "
        cmd += "--very-sensitive "
        cmd += "-x  {params.genome} "
        cmd += " {params.filetype} "

        if params.paired:
             cmd += '-1 {input[0]} -2 {input[1]} '
        else:
             cmd += '-U {input} '
        cmd += " 2> {log} | samtools view -bS - > {output} "
        shell(cmd)

INPUT_ALL.append(expand(rules.read_mapping.output, sample=samples.Name.tolist(), reference=config['reference']))

rule remove_chromosomes:
    input: join(OUT_DIR, '{sample}', '{reference}', 'mapping', 'sample.bam')
    output: join(OUT_DIR, '{sample}', '{reference}', 'mapping', 'sample.cleanchrom.bam'), \
            join(OUT_DIR, '{sample}', '{reference}', 'report', 'summary_removed_chroms.tsv')
    params: chroms = lambda wc: config['reference'][wc.reference]['removechroms']
    resources:
        mem_mb=1000
    run:
        remove_chroms(input[0], output[0], params.chroms, output[1])

INPUT_ALL.append(expand(rules.remove_chromosomes.output, sample=samples.Name.tolist(), reference=config['reference']))


rule make_genome_size_table:
    input: join(OUT_DIR, '{sample}', '{reference}', 'mapping', 'sample.cleanchrom.bam')
    output: join(OUT_DIR, '{sample}', '{reference}', 'mapping', '{reference}.genome')
    resources:
        mem_mb=1000
    run:
        make_genome_size(input[0], output[0])

INPUT_ALL.append(expand(rules.make_genome_size_table.output, reference=config['reference'], sample=samples.Name.tolist()))
