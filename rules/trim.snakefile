from os.path import splitext

if config['trim_reads'] == 'flexbar':
    rule adapter_trimming_flexbar:
        """Trim adapters using flexbar"""
        input:
            reads=get_trim_inputs
        params:
            target=TRIM_PATTERN,
            paired=is_paired
        output: TRIM_PATTERN + '_1.fastq', TRIM_PATTERN + '_2.fastq', TRIM_PATTERN + '.fastq'
        threads: 40
        log: join(LOG_DIR, 'flexbar_{sample}.log')
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
    
elif config['trim_reads'] == 'trim_galore':
    rule adapter_trimming_trimgalore:
        """Trim adapters using flexbar"""
        input:
            reads=get_trim_inputs
        params:
            target=TRIM_PATTERN,
            paired=is_paired,
            reads= get_trimgalore_output
        output: TRIM_PATTERN + '_1.fastq.gz', \
                TRIM_PATTERN + '_2.fastq.gz'
                #TRIM_PATTERN + '.fastq', directory(TRIM_PATTERN)
        threads: 40
        log: join(LOG_DIR, 'trimgalore_{sample}.log')
        run:
            cmd = 'mkdir -p {params.target} && trim_galore '
            if params.paired:
                cmd += " --paired {input.reads[0]} {input.reads[1]} "
            else:
                cmd += " {input.reads[0]} "
    
            cmd += " -o {params.target} --gzip"
            if 'adapters' in config:
                cmd += ' -a {}'.format(config['adapters'])
            cmd += " >> {log} 2>&1 && "
    
            # rename the trim galore output.
            if params.paired:
                cmd += " mv {params.target}/{params.reads[0]} {params.target}_1.fastq.gz; "
                cmd += " mv {params.target}/{params.reads[1]} {params.target}_2.fastq.gz; "
            else:
                cmd += " mv {params.target}/{params.reads[0]} {params.target}_1.fastq.gz; "

            # create fake output files
            # this allows us to use a single rule for 
            # single or paired end.
            if params.paired:
                cmd += " ln -s {params.target}_1.fastq.gz {params.target}.fastq.gz; "
            else:
                cmd += " ln -s {params.target}.fastq.gz {params.target}_1.fastq.gz; "
                cmd += " ln -s {params.target}.fastq.gz {params.target}_2.fastq.gz; "
            shell(cmd)
