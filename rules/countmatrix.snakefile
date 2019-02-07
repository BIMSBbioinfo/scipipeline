import traceback
import sys

# ------------------------- #
# Create counting bed files

rule make_counting_bins:
    input:
        bams = join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam"),
        bai = join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam.bai")
    output:
        bins = join(OUT_DIR, '{sample}', '{reference}', 'countmatrix', 'genomebins_binsize{binsize}.minmapq{minmapq}.mincount{mincounts}.bed')
    wildcard_constraints:
       minmapq='\d+', mincounts='\d+'
    resources:
        mem_mb=5000
    run:
        make_beds_for_intervalsize(input.bams, int(wildcards.binsize), output.bins)


# ------------------------- #
# Barcode-table

rule make_barcode_table:
    input:
        bams = join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.bam"),
    output:
        barcodes = join(OUT_DIR, '{sample}', '{reference}', 'countmatrix', 'sample_barcodes.tsv')
    resources:
        mem_mb=1000
    run:
        make_barcode_table(input.bams, output.barcodes)


# ------------------------- #
# Count reads in bins

rule counting_reads_in_bins:
    """Counting reads per barcode"""
    input:
        bams = join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam"),
        bai = join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam.bai"),
        bins = join(OUT_DIR, '{sample}', '{reference}', 'countmatrix', 'genomebins_binsize{binsize}.minmapq{minmapq}.mincount{mincounts}.bed')
    output:
        countmatrix = join(OUT_DIR, '{sample}', '{reference}', 'countmatrix', 'genomebins_binsize{binsize}.minmapq{minmapq}.mincount{mincounts}.tab'),
        cellsum = join(OUT_DIR, '{sample}', '{reference}', 'countmatrix', 'genomebins_binsize{binsize}.minmapq{minmapq}.mincount{mincounts}.tab.counts')
    params:
        count_both_ends=config['count_both_ends']
    wildcard_constraints:
       minmapq='\d+', mincounts='\d+'
    log: join(LOG_DIR, 'genomebins_countmatrix_{sample}_{reference}.binsize{binsize}.minmap{minmapq}.mincount{mincounts}.log')
    resources:
        mem_mb=10000
    run:
        try:
            sparse_count_reads_in_regions(input.bams, input.bins, output.countmatrix, flank=0, log=log[0], count_both_ends=params.count_both_ends)
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print("custom code exception:")
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                              limit=2, file=sys.stdout)


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
        bams = join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam"),
        bai = join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam.bai"),
        regions = join(OUT_DIR, '{sample}', "{reference}", "peaks", "sample.minmapq{minmapq}.mincount{mincounts}.flank{flank}_summits.bed")
    output: join(OUT_DIR, '{sample}', '{reference}', 'countmatrix', 'peakcounts_flank{flank}.minmapq{minmapq}.mincount{mincounts}.tab'),
            join(OUT_DIR, '{sample}', '{reference}', 'countmatrix', 'peakcounts_flank{flank}.minmapq{minmapq}.mincount{mincounts}.tab.counts')
    wildcard_constraints:
       minmapq='\d+', mincounts='\d+', flank='\d+'
    params:
        count_both_ends=config['count_both_ends']
    resources:
        mem_mb=10000
    log: join(LOG_DIR, 'peak_countmatrix_{sample}_{reference}.flank{flank}.minmap{minmapq}.mincount{mincounts}.log')
    run: 
        try:
            sparse_count_reads_in_regions(input.bams, input.regions, output[0], flank=0, log=log[0], count_both_ends=params.count_both_ends)
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print("custom code exception:")
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                              limit=2, file=sys.stdout)



INPUT_ALL.append(expand(rules.counting_reads_in_peaks.output,
                        reference=config['reference'],
                        sample=samples.Name.tolist(),
                        flank=config['peak_flank'],
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))


rule counting_reads_in_extra_annotation:
    """Counting reads per barcode"""
    input:
        bams = join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam"),
        bai = join(OUT_DIR, '{sample}', "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam.bai"),
        bins = lambda wc: config['reference'][wc.reference]['annotation'][wc.annotation]
    output:
        countmatrix = join(OUT_DIR, '{sample}', '{reference}', 'countmatrix', 'extra_{annotation}.minmapq{minmapq}.mincount{mincounts}.tab'),
        cellsum = join(OUT_DIR, '{sample}', '{reference}', 'countmatrix', 'extra_{annotation}.minmapq{minmapq}.mincount{mincounts}.tab.counts')
    wildcard_constraints:
       minmapq='\d+', mincounts='\d+'
    params:
        count_both_ends=config['count_both_ends']
    log: join(LOG_DIR, 'extra_countmatrix_{annotation}_{sample}_{reference}.minmap{minmapq}.mincount{mincounts}.log')
    resources:
        mem_mb=10000
    run:
        try:
            sparse_count_reads_in_regions(input.bams, input.bins, output.countmatrix, flank=0, log=log[0], count_both_ends=params.count_both_ends)
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print("custom code exception:")
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                              limit=2, file=sys.stdout)

for ref in config['reference']:
  if 'annotation' in config['reference'][ref]:
    print(ref)
    print(config['reference'][ref]['annotation'])
    INPUT_ALL.append(expand(rules.counting_reads_in_extra_annotation.output,
                            reference=[ref],
                            sample=samples.Name.tolist(),
                            minmapq=config['min_mapq'],
                            mincounts=config['min_counts_per_barcode'],
                            annotation=config['reference'][ref]['annotation']))
    
