from os.path import join

from utils.report_plots import plot_barcode_frequencies
from utils.report_plots import plot_fragment_size
from utils.count_matrix import get_barcode_frequency_genomewide


rule determine_barcode_frequencies:
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam")
    output: join(OUT_DIR, "{reference}", "report", "barcode_frequencies.{sample}.tab")
    run:
        get_barcode_frequency_genomewide(input[0], output[0])


# ------------------------ #
#
rule plot_barcode_freqs:
    """Plot barcode frequency"""
    input: join(OUT_DIR, "{reference}", "report", "barcode_frequencies.{sample}.tab")
    output: report(join(OUT_DIR, "{reference}", "report", "barcode_frequencies.{sample}.svg"))
    run:
        plot_barcode_frequencies(input[0], output[0])

INPUT_ALL.append(expand(rules.plot_barcode_freqs.output, reference=config['reference'], sample=samples.Name.tolist()))

rule plot_fragment_size_dist:
    """Plot fragment size distribution"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam")
    output: report(join(OUT_DIR, "{reference}", "{sample}.fragmentsize.svg"))
    run:
        plot_fragment_size(input[0], output[0])

INPUT_ALL.append(expand(rules.plot_fragment_size_dist.output, reference=config['reference'], sample=samples.Name.tolist()))


rule make_multiqc_report:
    input: expand(join(OUT_DIR, "{reference}", "macs2", "{sample}_summits.bed"), \
                  reference=config['reference'], sample=samples.Name.tolist()), \
           expand(join(OUT_DIR, 'fastqc', '{sample}'),  \
                  reference=config['reference'], sample=samples.Name.tolist())
    output: join(OUT_DIR, 'multiqc_report.html'), join(OUT_DIR, 'multiqc_data')
    params: searchdir=OUT_DIR
    shell:
        "multiqc {params.searchdir}"

INPUT_ALL.append(rules.make_multiqc_report.output)

