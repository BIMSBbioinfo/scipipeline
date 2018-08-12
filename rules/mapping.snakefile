from utils.cleanup_alignments import remove_chroms

def _bowtie_input_type_read(wildcards):
    filename = samples[samples.Name == wildcards.sample].read1.tolist()[0]
    return bowtie_input_filetype_option(filename)

rule read_mapping:
    "Maps reads against reference genome"
    input: get_mapping_inputs
    output: join(OUT_DIR, '{reference}', '{sample}.bam')
    wildcard_constraints: sample='\w+'
    params:
        genome=lambda wildcards: config['reference'][wildcards.reference]['bowtie2index'],
        paired=is_paired,
        filetype=_bowtie_input_type_read
    threads: 20
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
    input: join(OUT_DIR, '{reference}', '{sample}.bam')
    output: join(OUT_DIR, '{reference}', '{sample}.cleanchrom.bam')
    params: chroms = lambda wc: config['reference'][wc.reference]['removechroms']
    run:
        remove_chroms(input[0], output[0], params.chroms)

INPUT_ALL.append(expand(rules.remove_chromosomes.output, sample=samples.Name.tolist(), reference=config['reference']))
