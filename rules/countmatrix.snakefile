

# ------------------------- #
# Count reads in bins

rule counting_reads_in_bins:
    """Counting reads per barcode"""
    input:
        bams = join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam"),
        bai = join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam.bai")
    output:
        bins = join(OUT_DIR, '{reference}', 'countmatrix', '{sample}_bin_{binsize}.bed'),
        countmatrix = join(OUT_DIR, '{reference}', 'countmatrix', '{sample}_bin_{binsize}.tab')
    run:
        make_beds_for_intervalsize(input.bams, int(wildcards.binsize), output.bins)
        sparse_count_reads_in_regions(input.bams, output.bins, output.countmatrix, flank=0)


INPUT_ALL.append(expand(rules.counting_reads_in_bins.output,
                        reference=config['reference'],
                        binsize=config['binsize'],
                        sample=config['samples'].keys()))

# ------------------------- #
# Count reads in regions

rule counting_reads_in_peaks:
    """Counting reads per barcode"""
    input:
        bams = join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam"),
        bai = join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam.bai"),
        regions = join(OUT_DIR, "{reference}", "macs2", "{sample}_summits.bed")
    output: join(OUT_DIR, '{reference}', 'countmatrix', '{sample}_peak_counts_flank_{flank}.h5')
    run: sparse_count_reads_in_regions(input.bams, input.regions, \
         output[0], flank=int(wildcards.flank))


INPUT_ALL.append(expand(rules.counting_reads_in_peaks.output,
                        reference=config['reference'],
                        sample=config['samples'].keys(),
                        flank=config['peak_flank']))