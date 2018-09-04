import numpy as np
import pysam
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns

def plot_barcode_frequencies(tab_file, plotname):
    """ Plot barcode frequency distribution"""
    x=pd.read_csv(tab_file, sep='\t')
    f = plt.figure()
    ax = sns.distplot(x.counts.apply(np.log10))
    ax.set_xlabel('Log10(# Fragments)')
    ax.set_ylabel('Frequency')
    ax.set_title('Barcode frequency')
    f.savefig(plotname, dpi=f.dpi)



def plot_fragment_size(bamin, plotname):
    """Plot fragment size distribution"""
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
    """2D plot of percentage of reads in peaks vs. reads per barcode"""
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


def scatter_log_frequencies_per_species(tables, labels, plotname):
    t1 = pd.read_csv(tables[0], sep='\t')
    t2 = pd.read_csv(tables[1], sep='\t')
    for df in [t1, t2]:
        df['I1'] = df['file'].apply(lambda x: int(x.split('_')[0].split('-')[1]))

    #joined = pd.concat([t1,t2], axis=1, join='inner', on='file')
    joined = pd.merge(t1, t2, how='inner', on='file')
    joined['color'] = joined.I1_x.apply(lambda x: 'red' if x<5 else 'blue')
    joined = joined[['deduplicated_x', 'deduplicated_y', 'color']]
    joined.columns = labels + ['color']
    f, ax = plt.subplots()
    ax = joined.plot.scatter(labels[0], labels[1], ax=ax, loglog=True, alpha=.2, color=joined.color)
    ax.set_xlabel('Reads per {} barcode'.format(labels[0]))
    ax.set_ylabel('Reads per {} barcode'.format(labels[1]))
    #ax.set_xlabel(
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


def species_specificity(aln_file1, aln_file2, output, labels):
    aln1 = AlignmentFile(aln_file1, 'r').fetch(until_eof=True)
    aln2 = AlignmentFile(aln_file2, 'r').fetch(until_eof=True)

    cnt = np.zeros((2,2))
    try:
        while 1:
            m1 = not next(aln1).is_unmapped
            m2 = not next(aln2).is_unmapped
            cnt[0 if m1 else 1, 0 if m2 else 1] += 1
    except StopIteration:
        pass
    finally:
        print(cnt)
        f = plt.figure()
        sns.heatmap(cnt/cnt.sum(), annot=True, cmap='YlGnBu',
                    xticklabels=['mapped', 'unmapped'], yticklabels=['mapped', 'unmapped'])
        plt.ylabel(labels[0])
        plt.xlabel(labels[1])
        f.savefig(output + '.png', dpi=f.dpi)
        np.savetxt(output + '.tab', cnt, delimiter='\t')

