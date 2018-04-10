from os.path import join

from utils.report_plots import plot_barcode_frequencies
from utils.report_plots import species_specificity
from utils.report_plots import scatter_log_frequencies_per_species
from utils.report_plots import scatter_frequencies_per_species_colored
from utils.report_plots import density_frequencies_per_species_colored

# Reference genome and mapping index
BAM_SORTED = join(OUT_DIR, '{reference}', 'atac_sorted.bam')


# ------------------------ #
#
rule plot_barcode_freqs:
    """Plot barcode frequency"""
    input: join(OUT_DIR, "report", "barcode_frequencies.{reference}.tab")
    output: join(OUT_DIR, "report", "barcode_frequencies_{reference}.png")
    run:
        plot_barcode_frequencies(input[0], output[0])

INPUT_ALL.append(expand(rules.plot_barcode_freqs.output, reference=config['reference']))


# ------------------------ #
#
rule read_cooccurrence_by_species:
    """How frequently do reads map uniquely to one species or to both?"""
    input: expand(BAM_SORTED, reference=config['reference'])
    output: expand(join(OUT_DIR, "report", "read_cooccurrence.{suffix}"), suffix=['tab', 'png'])
    run:
        species_specificity(input[0], input[1], output[0].split('.')[0], expand('{reference}', reference=config['reference']))

INPUT_ALL.append(rules.read_cooccurrence_by_species.output)

# ------------------------ #
#
rule scatterplot_barcode_logfreqs:
    """Plot barcode frequency"""
    input: 
        tables = expand(join(OUT_DIR, "report", 
           "barcode_frequencies.{reference}.tab"), reference=config['reference']),
    params:
        names = expand('{reference}', reference=config['reference'])

    output: join(OUT_DIR, "report", "barcode_log_frequencies_scatter.png")
    run:
        scatter_log_frequencies_per_species(input.tables, params.names, output[0])

INPUT_ALL.append(rules.scatterplot_barcode_logfreqs.output)

# ------------------------ #
#
rule scatterplot_barcode_freqs:
    """Plot barcode frequency"""
    input: 
        tables = expand(join(OUT_DIR, "report", 
           "barcode_frequencies.{reference}.tab"), reference=config['reference']),
    params:
        names = expand('{reference}', reference=config['reference'])

    output: join(OUT_DIR, "report", "barcode_frequencies_scatter.png")
    run:
        scatter_frequencies_per_species_colored(input.tables, params.names, output[0])

INPUT_ALL.append(rules.scatterplot_barcode_freqs.output)

# ------------------------ #
#
rule density_barcode_freqs:
    """Plot barcode frequency"""
    input: 
        tables = expand(join(OUT_DIR, "report", 
           "barcode_frequencies.{reference}.tab"), reference=config['reference']),
    params:
        names = expand('{reference}', reference=config['reference'])

    output: join(OUT_DIR, "report", "barcode_frequencies_density.png")
    run:
        density_frequencies_per_species_colored(input.tables, params.names, output[0])

INPUT_ALL.append(rules.density_barcode_freqs.output)

