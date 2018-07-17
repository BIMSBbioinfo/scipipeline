

# ------------------------- #
# report barcode frequencies

rule peak_calling_on_aggregate:
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam")
    output: join(OUT_DIR, "{reference}", "macs2", "{sample}_peaks.narrowPeak"), join(OUT_DIR, "{reference}", "macs2", "{sample}_summits.bed")
    params: name='{sample}',
            outdir = join(OUT_DIR, "{reference}", "macs2"),
            foption = lambda wc: 'BAMPE' if is_paired(wc) else 'BAM',
            gsize = lambda wc: config['reference'][wc.reference]['macs_gsize']
    log: join(LOG_DIR, 'macs2_{sample}.log')
    shell:
      " macs2 callpeak --name {params.name} -t {input} -f " +
      "{params.foption}" +
      " --nomodel --keep-dup all --outdir {params.outdir} --call-summits --gsize {params.gsize} 2> {log} "

INPUT_ALL.append(expand(rules.peak_calling_on_aggregate.output, reference=config['reference'], sample=samples.Name.tolist()))
