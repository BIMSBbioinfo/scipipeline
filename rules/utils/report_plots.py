import numpy as np
import pysam
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.misc import comb
import seaborn as sns
sns.set(style='white')

def plot_barplot_summary_statistics(tab_file, plotname, title):
    """Plot summary statistics of filtering step."""
    df=pd.read_table(tab_file)
    f, ax = plt.subplots()
    ax = df.plot.barh(x='Readgroup', y='counts', ax=ax)
    ax.set_title(title)
    f.savefig(plotname, dpi=f.dpi)


def plot_barcode_frequencies(tab_file, plotname):
    """ Plot barcode frequency distribution.

    This function takes a tsv-file as input.
    The Tsv file consists of two columns, barcode names and read counts

    barcodes    counts
    bc1         502
    bc2         1304
    ...

    The output will be a figure showing the distribution of
    log10 transformed counts that is saved under plotname.

    Parameters
    ----------
    tab_file : str
        TSV table containing the barcodes with the associated counts.
    plotname : str
        Output filename of the figure.
    """
    x=pd.read_csv(tab_file, sep='\t')
    f = plt.figure()
    ax = sns.distplot(x.counts.apply(lambda x: np.log10(x + 1.)))
    ax.set_xlabel('Log10(# Fragments)')
    ax.set_ylabel('Frequency')
    ax.set_title('Barcode frequency')
    f.savefig(plotname, dpi=f.dpi)



def plot_fragment_size(bamin, plotname):
    """Plot fragment size distribution.

    This plot illustrates the fragment size distribution
    across all cells. It should show a periodic pattern
    corresponding to nucleosome free and nucleosome spanning
    read pairs.
    """
    handle = pysam.AlignmentFile(bamin, 'r')
    fragmentsize_dist = np.zeros((2000,))

    for aln in handle:
       if not aln.is_unmapped and aln.is_read1:
           tl = min(abs(aln.tlen), 1999)
           fragmentsize_dist[tl] += 1

    handle.close()

    # make a plot
    f = plt.figure()
    plt.plot(list(range(2000)), fragmentsize_dist)
    plt.title('Fragment size histogram')
    plt.xlabel('Fragment size')
    plt.ylabel('Frequency')
    f.savefig(plotname, dpi=f.dpi)


def plot_barcode_frequency_by_peak_percentage(barcode_frequency,
                                              peak_frequency, plotname):
    """2D plot of percentage of reads in peaks vs. reads per barcode.

    This plot shows distribution of fragments per cell relative
    to the percentage of reads falling into peaks.

    Buenrostro et al. 2018 have used this information to filter out
    low quality cells.
    """
    y = pd.read_csv(peak_frequency, sep='\t')
    y = y.groupby('cell').sum()['count']

    x = pd.read_csv(barcode_frequency, sep='\t')
    x = x.groupby('cell').sum()['count']
    y = y/x * 100.

    x = x.apply(np.log10)

    ax = sns.jointplot(x, y, kind='hex', ylim=[0, 100])
    ax = ax.set_axis_labels('Log10(# Fragments per cell)',
                            'Fragments in peaks [ % ]')
    ax = ax.savefig(plotname)


def barcode_collision_scatter_plot(tables, labels, plotname, logplot=True):
    """ Plots the barcode frequencies per species.

    This function gets a list of barcode frequency tables from different species
    and plots the frequencies of one species against another one.
    This figure will indicate issues barcode collisions.

    Parameters
    ----------
    tables : list(str)
        List of Tsv files containing the barcode counts. Each file corresponds
        to a species with barcodes from the same barcode universe.
    labels : list(str)
        List of species labels.
    plotname : str
        Figure name. If more than 2 species are used, the figure will contain
        subplots for each pair-wise comparison.
    logplot : boolean
        Whether to show log transformed counts in the scatter plot or raw counts.
        Default: True.
    """
    nfigures = int(comb(len(tables), 2))

    f, axes = plt.subplots(int(np.ceil(nfigures / 2)), 2)

    #print(len(axes), len(axes[0]))
    print(labels, len(labels))
    print(tables, len(tables))

    for xdim in range(len(tables) - 1):
        for ydim in range(xdim + 1, len(tables)):

            t1 = pd.read_csv(tables[xdim], sep='\t')
            t2 = pd.read_csv(tables[ydim], sep='\t')

            joined = pd.merge(t1, t2, how='inner', on='barcodes')
            print(joined.head())

            axes[xdim, ydim - xdim - 1] = joined.plot.scatter('counts_x',
                                                   'counts_y',
                                                   ax=axes[xdim, ydim - xdim - 1],
                                                   loglog=logplot,
                                                   alpha=.2)
            labtext = "{} for ".format("Log10(#Fragments)" if logplot else "#Fragments")

            axes[xdim, ydim - xdim - 1].set_xlabel(labtext + labels[xdim])
            axes[xdim, ydim - xdim - 1].set_ylabel(labtext + labels[ydim])

    f.savefig(plotname, dpi=f.dpi)


def scatter_frequencies_per_species_colored(tables, labels, plotname):
    t1 = pd.read_csv(tables[0], sep='\t')
    t2 = pd.read_csv(tables[1], sep='\t')
    for df in [t1, t2]:
        df['I1'] = df['file'].apply(lambda x: int(x.split('_')[0].split('-')[1]))

    #joined = pd.concat([t1,t2], axis=1, join='inner', on='file')
    joined = pd.merge(t1, t2, how='inner', on='file')
    joined['color'] = joined.I1_x.apply(lambda x: 'red' if x<5 else 'blue')
    joined = joined[['deduplicated_x', 'deduplicated_y', 'color']]
    joined.columns = labels + ['color']
    #joined = joined.apply(np.log10)
    f, ax = plt.subplots()
    ax = joined.plot.scatter(labels[0], labels[1], ax=ax, alpha=.2, color=joined.color)
    ax.set_xlabel('Reads per {} barcode'.format(labels[0]))
    ax.set_ylabel('Reads per {} barcode'.format(labels[1]))
    ax.set_xlim(0, 70000)
    ax.set_ylim(0, 70000)
    #ax.set_xlabel(
    f.savefig(plotname, dpi=f.dpi)


def density_frequencies_per_species_colored(tables, labels, plotname):
    t1 = pd.read_csv(tables[0], sep='\t')
    t2 = pd.read_csv(tables[1], sep='\t')
    for df in [t1, t2]:
        df['I1'] = df['file'].apply(lambda x: int(x.split('_')[0].split('-')[1]))

    joined = pd.merge(t1, t2, how='inner', on='file')
    joined['color'] = joined.I1_x.apply(lambda x: 'red' if x<5 else 'blue')
    joined = joined[['deduplicated_x', 'deduplicated_y', 'color']]
    joined.columns = labels + ['color']

    f, ax = plt.subplots()
    ax.set_aspect('equal')

    fish_barcodes = joined.query("color == 'red'")
    fly_barcodes = joined.query("color == 'blue'")
    ax = sns.kdeplot(fish_barcodes[labels[0]],
                     fish_barcodes[labels[1]],
                     cmap="Reds", shade=True, shade_lowest=False, alpha=.3)
    ax = sns.kdeplot(fly_barcodes[labels[0]],
                     fly_barcodes[labels[1]],
                     cmap="Blues", shade=True, shade_lowest=False, alpha=.3)
    ax.set_xlabel('Reads per {} barcode'.format(labels[0]))
    ax.set_ylabel('Reads per {} barcode'.format(labels[1]))
    ax.set_xlim(0, 70000)
    ax.set_ylim(0, 70000)
    red = sns.color_palette("Reds")[-2]
    blue = sns.color_palette("Blues")[-2]
    f.savefig(plotname, dpi=f.dpi)


def cross_species_mapping_reads(bamfiles, plotname, labels):
    """This function determines the cross-mapping reads.

    In order to determine the success of a combinatorial indexing
    experiment, usually one mixes reads obtained from two different
    species. This allows to investigate and detect issues with
    barcode collisions.
    However, barcode may also occur in two species if they in fact
    map to both species simultaneously through cross mapability.
    This function produces a heatmap that visualizes how many reads
    map to one or more than one species.

    Parameters
    ----------
    bamfiles : list(str)
        List of bam files
    plotname : str
        Output location of figure.
    labels : list(str)
        List of labels included in the figure
    """

    readers = [AlignmentFile(file_, 'r').fetch(until_eof=True) for file_ in bamfiles]

    cnt = np.zeros((len(readers), len(readers)))
    try:
        while 1:
            vals = [not next(aln).is_unmapped for aln in readers]
            for i, vali in enumerate(vals[:-1]):
                if vali == 0:
                    continue
                for j, valj in enumerate(vals[1:]):
                    if valj == 0:
                        continue
                    # if vali and valj are both 1, set the value in the
                    # confusion matrix
                    cnt[i, j] += 1
    except StopIteration:
        pass
    finally:
        f = plt.figure()
        sns.heatmap((cnt/cnt.sum(axis=0)).T, annot=True, cmap='YlGnBu',
                    xticklabels=labels, yticklabels=labels)
        plt.ylabel("")
        plt.xlabel("(Row-normal.) Cross mappability")
        f.savefig(plotname, dpi=f.dpi)
