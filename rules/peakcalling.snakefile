

# ------------------------- #
# report barcode frequencies

rule peak_calling_on_aggregate:
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    output: join(OUT_DIR, "{reference}", "macs2", "{sample}.minmapq{minmapq}.mincount{mincounts}_peaks.narrowPeak"), join(OUT_DIR, "{reference}", "macs2", "{sample}.minmapq{minmapq}.mincount{mincounts}_summits.bed")
    params: name='{sample}.minmapq{minmapq}.mincount{mincounts}',
            outdir = join(OUT_DIR, "{reference}", "macs2"),
            foption = lambda wc: 'BAMPE' if is_paired(wc) else 'BAM',
            gsize = lambda wc: config['reference'][wc.reference]['macs_gsize']
    log: join(LOG_DIR, 'macs2_{sample}_{reference}_minmapq{minmapq}_mincount{mincounts}.log')
    shell:
      " macs2 callpeak --name {params.name} -t {input} -f " +
      "{params.foption}" +
      " --nomodel --keep-dup all --outdir {params.outdir} --call-summits --gsize {params.gsize} 2> {log} "

INPUT_ALL.append(expand(rules.peak_calling_on_aggregate.output, 
                        reference=config['reference'], 
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))
