
rule read_mapping:
    "Maps reads against reference genome"
    input: get_mapping_inputs
    output: temp(join(OUT_DIR, '{reference}', '{sample}.bam'))
    params:
        genome=lambda wildcards: config['reference'][wildcards.reference]['bowtie2index'],
        paired=is_paired,
        filetype = lambda wildcards: bowtie_input_filetype_option(config['samples'][wildcards.sample]['read1'])
    threads: 20
    log: join(LOG_DIR, '{sample}_{reference}_bowtie2.log')
    run:
        cmd = 'bowtie2'
        cmd = "bowtie2 -p {threads} -X 2000 --no-mixed --no-discordant "
        cmd += "-x  {params.genome} "
        cmd += " {params.filetype} "

        if params.paired:
             cmd += '-1 {input[0]} -2 {input[1]} '
        else:
             cmd += '-U {input} '
        cmd += " 2> {log} | samtools view -bS - > {output} "
        shell(cmd)

INPUT_ALL.append(expand(rules.read_mapping.output, sample=config['samples'].keys(), reference=config['reference']))
