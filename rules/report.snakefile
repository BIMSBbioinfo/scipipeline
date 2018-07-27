from os.path import join

from utils.report_plots import plot_barcode_frequencies
from utils.report_plots import plot_fragment_size


# ------------------------ #
#
rule plot_barcode_freqs:
    """Plot barcode frequency"""
    input: join(OUT_DIR, "report", "barcode_frequencies.{reference}.tab")
    output: report(join(OUT_DIR, "report", "barcode_frequencies_{reference}.svg"))
    run:
        plot_barcode_frequencies(input[0], output[0])

INPUT_ALL.append(expand(rules.plot_barcode_freqs.output, reference=config['reference']))

rule fragment_size_distribution:
    """Plot fragment size distribution"""
    input: join(OUT_DIR, "{reference}", "{sample}.barcoded.dedup.bam")
    output: report(join(OUT_DIR, "{reference}", "{sample}.fragmentsize.svg"))
    run
        plot_fragment_size(input[0], output[0])

INPUT_ALL.append(expand(rules.plot_fragment_size.output, reference=config['reference'], sample=samples.Name.tolist()))
