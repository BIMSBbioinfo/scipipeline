

# ------------------------- #
# report barcode frequencies

rule peak_calling_on_aggregate:
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    output: join(OUT_DIR, "{reference}", "macs2", "{sample}.minmapq{minmapq}.mincount{mincounts}_peaks.narrowPeak"), \
            join(OUT_DIR, "{reference}", "macs2", "{sample}.minmapq{minmapq}.mincount{mincounts}_summits.bed")
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

rule merge_overlapping_peaks:
    """ Peak summits extended by flank, merged, sorted and trimmed to same size"""
    input: join(OUT_DIR, "{reference}", "macs2", "{sample}.minmapq{minmapq}.mincount{mincounts}_summits.bed")
    output: join(OUT_DIR, "{reference}", "macs2", "{sample}.minmapq{minmapq}.mincount{mincounts}.flank{flk}_summits.bed")
    params:
       genome = lambda wc: config['reference'][wc.reference]['genomesize'],
       flks = lambda wc: wc.flk
    shell:
      """
      bedtools flank -i {input} -g {params.genome} -l {params.flks} -r {params.flks} | bedtools sort -i stdin | bedtools merge -i stdin | awk -F '{{ OFS="\\t"; mid=($2+$3)/2; print $1, mid-{params.flks},mid+{params.flks}}}' > {output}"
     """


INPUT_ALL.append(expand(rules.merge_overlapping_peaks.output, 
                        reference=config['reference'], 
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        flk=config['peak_flank'],
                        mincounts=config['min_counts_per_barcode']))

