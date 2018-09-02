from os.path import join

from utils.report_plots import plot_barcode_frequencies
from utils.report_plots import plot_fragment_size
from utils.count_matrix import get_barcode_frequency_genomewide


rule determine_barcode_frequencies:
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    output: join(OUT_DIR, "{reference}", "report", "barcode_frequencies.{sample}.minmapq{minmapq}.mincount{mincounts}.tab")
    run:
        get_barcode_frequency_genomewide(input[0], output[0])


# ------------------------ #
#
rule plot_barcode_freqs:
    """Plot barcode frequency"""
    input: join(OUT_DIR, "{reference}", "report", "barcode_frequencies.{sample}.minmapq{minmapq}.mincount{mincounts}.tab")
    output: report(join(OUT_DIR, "{reference}", "report", "barcode_frequencies.{sample}.minmapq{minmapq}.mincount{mincounts}.svg"))
    run:
        plot_barcode_frequencies(input[0], output[0])

INPUT_ALL.append(expand(rules.plot_barcode_freqs.output, reference=config['reference'], sample=samples.Name.tolist(), 
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))

rule plot_fragment_size_dist:
    """Plot fragment size distribution"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    output: report(join(OUT_DIR, "{reference}", "report", "{sample}.fragmentsize_minmapq{minmapq}_mincount{mincounts}.svg"))
    run:
        plot_fragment_size(input[0], output[0])

INPUT_ALL.append(expand(rules.plot_fragment_size_dist.output, 
                        reference=config['reference'], 
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))


rule make_multiqc_report:
    input: expand(join(OUT_DIR, "{reference}", "macs2", \
                       "{sample}.minmapq{minmapq}.mincount{mincounts}_summits.bed"), \
                  reference=config['reference'], \
                  sample=samples.Name.tolist(), \
                  minmapq=config['min_mapq'], \
                  mincounts=config['min_counts_per_barcode'])
    output: join(OUT_DIR, 'multiqc_report.html'), directory(join(OUT_DIR, 'multiqc_data'))
    params: searchdir=OUT_DIR
    shell:
        "multiqc --outdir {params.searchdir} {params.searchdir}"

INPUT_ALL.append(rules.make_multiqc_report.output)

