

# ------------------------- #
# report barcode frequencies


rule sort_for_peak_calling_on_aggregate:
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    output: temp(join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}_namesort.bam"))
    resources:
      mem_mb=30000
    shell:
        "samtools sort -n {input} -o {output}"


rule peak_calling_on_aggregate:
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}_namesort.bam")
    output: join(OUT_DIR, "{sample}", "{reference}", "peaks", "sample.minmapq{minmapq}.mincount{mincounts}_peaks.narrowPeak")
    resources:
      mem_mb=3000
    log: join(LOG_DIR, 'macs2_{sample}_{reference}_minmapq{minmapq}_mincount{mincounts}.log')
    shell:
      " Genrich -t {input} -o {output} -j -y -v -g 175 -p 0.05  2> {log} "

INPUT_ALL.append(expand(rules.peak_calling_on_aggregate.output, 
                        reference=config['reference'], 
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))

rule merge_overlapping_peaks:
    """ Peak summits extended by flank, merged, sorted and trimmed to same size"""
    input: 
       inbed = join(OUT_DIR, "{sample}", "{reference}", "peaks", "sample.minmapq{minmapq}.mincount{mincounts}_peaks.narrowPeak"), \
       genomesize = join(OUT_DIR, '{sample}', "{reference}", 'mapping', "{reference}.genome")
    output: join(OUT_DIR, "{sample}", "{reference}", "peaks", "sample.minmapq{minmapq}.mincount{mincounts}.flank{flk}_merged.bed")
    params:
       flks = lambda wc: wc.flk
    resources:
       mem_mb=2000
    shell:
      """
      bedtools slop -i {input.inbed} -g {input.genomesize} -l {params.flks} -r {params.flks} | bedtools sort -i stdin | bedtools merge -i stdin | awk '{{ OFS="\\t"; mid=int(($2+$3)/2); print $1, mid-{params.flks},mid+{params.flks}}}' > {output}
     """


INPUT_ALL.append(expand(rules.merge_overlapping_peaks.output, 
                        reference=config['reference'], 
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        flk=config['peak_flank'],
                        mincounts=config['min_counts_per_barcode']))

