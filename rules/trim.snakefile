

rule adapter_trimming:
    """Trim adapters using flexbar"""
    input:
        reads=get_trim_inputs
    params:
        target=TRIM_PATTERN,
        paired=is_paired
    output: TRIM_PATTERN + '_1.fastq', TRIM_PATTERN + '_2.fastq', TRIM_PATTERN + '.fastq'
    threads: 40
    log: join(LOG_DIR, 'flexbar_{sample}.log')
    message:"""
        Trimming fastq files:
            input: {input}
            output: {output}
    """
    run:
        cmd = 'flexbar '
        cmd += "-r {input.reads[0]} "
        if params.paired:
            cmd += " -p {input.reads[1]} "
        cmd += " -t {params.target}"
        if 'adapters' in config:
            cmd += ' -a {}'.format(config['adapters'])
        cmd += " -f i1.8 -u 10 -ae RIGHT -at 1.0 --threads {threads} "
        cmd += " --min-read-length 50  > {log} && "
        if params.paired:
            cmd += " ln -s {params.target}_1.fastq {params.target}.fastq; "
        else:
            # create fake output files for single-end data
            cmd += " ln -s {params.target}.fastq {params.target}_1.fastq; "
            cmd += " ln -s {params.target}.fastq {params.target}_2.fastq; "
        shell(cmd)
