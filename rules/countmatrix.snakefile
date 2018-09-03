

# ------------------------- #
# Create counting bed files

rule make_counting_bins:
    input:
        bams = join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam"),
        bai = join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam.bai")
    output:
        bins = join(OUT_DIR, '{reference}', 'countmatrix', 'genomebins_{sample}_binsize{binsize}.minmapq{minmapq}.mincount{mincounts}.bed')
    wildcard_constraints:
       minmapq='\d+', mincounts='\d+'
    run:
        make_beds_for_intervalsize(input.bams, int(wildcards.binsize), output.bins)


# ------------------------- #
# Barcode-table

rule make_barcode_table:
    input:
        bams = join(OUT_DIR, "{reference}", "{sample}.barcoded.bam"),
    output:
        barcodes = join(OUT_DIR, '{reference}', 'countmatrix', '{sample}_barcodes.tsv')
    run:
        make_barcode_table(input.bams, output.barcodes)


# ------------------------- #
# Count reads in bins

rule counting_reads_in_bins:
    """Counting reads per barcode"""
    input:
        bams = join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam"),
        bai = join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam.bai"),
        bins = join(OUT_DIR, '{reference}', 'countmatrix', 'genomebins_{sample}_binsize{binsize}.minmapq{minmapq}.mincount{mincounts}.bed')
    output:
        countmatrix = join(OUT_DIR, '{reference}', 'countmatrix', 'genomebins_{sample}_binsize{binsize}.minmapq{minmapq}.mincount{mincounts}.tab'),
        cellsum = join(OUT_DIR, '{reference}', 'countmatrix', 'genomebins_{sample}_binsize{binsize}.minmapq{minmapq}.mincount{mincounts}.tab.counts')
    wildcard_constraints:
       minmapq='\d+', mincounts='\d+'
    run:
        sparse_count_reads_in_regions(input.bams, input.bins, output.countmatrix, flank=0)


INPUT_ALL.append(expand(rules.counting_reads_in_bins.output,
                        reference=config['reference'],
                        binsize=config['binsize'],
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))

# ------------------------- #
# Count reads in regions

rule counting_reads_in_peaks:
    """Counting reads per barcode"""
    input:
        bams = join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam"),
        bai = join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam.bai"),
        regions = join(OUT_DIR, "{reference}", "macs2", "{sample}.minmapq{minmapq}.mincount{mincounts}.flank{flank}_summits.bed")
    output: join(OUT_DIR, '{reference}', 'countmatrix', 'peak_counts_{sample}_flank{flank}.minmapq{minmapq}.mincount{mincounts}.tab'),
            join(OUT_DIR, '{reference}', 'countmatrix', 'peak_counts_{sample}_flank{flank}.minmapq{minmapq}.mincount{mincounts}.tab.counts')
    wildcard_constraints:
       minmapq='\d+', mincounts='\d+', flank='\d+'
    run: sparse_count_reads_in_regions(input.bams, input.regions, output[0])


INPUT_ALL.append(expand(rules.counting_reads_in_peaks.output,
                        reference=config['reference'],
                        sample=samples.Name.tolist(),
                        flank=config['peak_flank'],
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))
