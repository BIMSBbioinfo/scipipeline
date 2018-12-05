from os.path import join

from utils.report_plots import plot_barcode_frequencies
from utils.report_plots import plot_barcode_frequency_by_peak_percentage
from utils.report_plots import plot_fragment_size
from utils.count_matrix import get_barcode_frequency_genomewide


# ------------------------ #
#
rule determine_barcode_frequencies:
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    output: join(OUT_DIR, "{sample}", "{reference}", "report", "barcode_frequencies.minmapq{minmapq}.mincount{mincounts}.tab")
    resources:
      mem_mb=1000
    run:
        get_barcode_frequency_genomewide(input[0], output[0])


# ------------------------ #
#
rule plot_barcode_freqs:
    """Plot barcode frequency"""
    input: join(OUT_DIR, "{sample}", "{reference}", "report", "barcode_frequencies.minmapq{minmapq}.mincount{mincounts}.tab")
    output: report(join(OUT_DIR, "{sample}", "{reference}", "report", "barcode_frequencies.minmapq{minmapq}.mincount{mincounts}.svg"), category="Barcode frequency")
    resources:
      mem_mb=1000
    run:
        plot_barcode_frequencies(input[0], output[0])

INPUT_ALL.append(expand(rules.plot_barcode_freqs.output, reference=config['reference'], sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))

# ------------------------ #
#
rule plot_barcode_freqs_by_percent_in_peaks:
    """Plot barcode frequency by percent in peaks"""
    input: x = expand(join(OUT_DIR, '{{sample}}', '{{reference}}', 'countmatrix', 'genomebins_binsize{binsize}.minmapq{{minmapq}}.mincount{{mincounts}}.tab'), binsize=config['binsize']), \
           y = join(OUT_DIR, "{sample}", '{reference}', 'countmatrix', 'peakcounts_flank{flank}.minmapq{minmapq}.mincount{mincounts}.tab')
    output: report(join(OUT_DIR, "{sample}", "{reference}", "report", "reads_in_peak_percentage.minmapq{minmapq}.mincount{mincounts}.flank{flank}.binsize{binsize}.svg"), category="Barcode frequency")
    resources:
      mem_mb=4000
    run:
        # input.x as a list containing different binsizes.
        # however, it is sufficient to just use one of them.
        plot_barcode_frequency_by_peak_percentage(input.x[0], input.y, output[0])

INPUT_ALL.append(expand(rules.plot_barcode_freqs_by_percent_in_peaks.output,
                        reference=config['reference'],
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode'],
                        binsize=config['binsize'],
                        flank=config['peak_flank']))

# ------------------------ #
#
rule plot_fragment_size_dist:
    """Plot fragment size distribution"""
    input: join(OUT_DIR, "{sample}", "{reference}", 'mapping', "sample.barcoded.minmapq{minmapq}.dedup.mincount{mincounts}.bam")
    output: report(join(OUT_DIR, "{sample}", "{reference}", "report", "fragmentsizedist_minmapq{minmapq}_mincount{mincounts}.svg"), category="Fragment size distribution")
    resources:
      mem_mb=1000
    run:
        plot_fragment_size(input[0], output[0])

INPUT_ALL.append(expand(rules.plot_fragment_size_dist.output,
                        reference=config['reference'],
                        sample=samples.Name.tolist(),
                        minmapq=config['min_mapq'],
                        mincounts=config['min_counts_per_barcode']))


# ------------------------ #
#
if len(config['reference']) > 1:
    # if more than one species is available, a figure is created
    # to check for barcode collisions.
    rule plot_barcode_collision_across_species:
        """Plot fragment size distribution"""
        input:
            tables = expand(join(OUT_DIR, "{{sample}}", \
                        "{reference}", "report", \
                        "barcode_frequencies.minmapq{{minmapq}}.mincount{{mincounts}}.tab"), \
                        reference=config['reference'])
        params:
            names = config['reference']
        output: report(join(OUT_DIR, "{sample}", "report", \
                       "barcode_collision_logplot_minmapq{minmapq}_mincount{mincounts}.svg"), \
                       category="Barcode collision"),
                report(join(OUT_DIR, "{sample}", "report", \
                               "barcode_collision_plot_minmapq{minmapq}_mincount{mincounts}.svg"), \
                               category="Barcode collision"),
        resources:
          mem_mb=1000
        run:
            barcode_collision_scatter_plot(input.tables, params.names, \
                                                output[0], logplot=True)
            barcode_collision_scatter_plot(input.tables, params.names, \
                                                output[1], logplot=False)

    INPUT_ALL.append(expand(rules.plot_barcode_collision_across_species.output,
                            sample=samples.Name.tolist(),
                            minmapq=config['min_mapq'],
                            mincounts=config['min_counts_per_barcode']))


# ------------------------ #
#
rule make_multiqc_report:
    input: expand(join(OUT_DIR, "{sample}", "{reference}", "peaks", \
                       "sample.minmapq{minmapq}.mincount{mincounts}_summits.bed"), \
                  reference=config['reference'], \
                  sample=samples.Name.tolist(), \
                  minmapq=config['min_mapq'], \
                  mincounts=config['min_counts_per_barcode'])
    output: join(OUT_DIR, 'multiqc_report.html'), directory(join(OUT_DIR, 'multiqc_data'))
    params: searchdir=OUT_DIR
    resources:
      mem_mb=10000
    shell:
        "multiqc -f --outdir {params.searchdir} {params.searchdir}"

INPUT_ALL.append(rules.make_multiqc_report.output)
