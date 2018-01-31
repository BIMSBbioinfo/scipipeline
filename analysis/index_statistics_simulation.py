"""Test index frequency.

The frequency with which indices occur in the dataset depends on
the cross-linking efficiency.
Ideally, if cross-linking works perfectly, all DNA fragments in a given
nucleus share the same barcode. In this case, the single-cell stage can
be derived.
However, if cross-linking is not working, the DNA fragments in a given
nucleus DO NOT share the same barcode. In this case, the cell state cannot
be reconstructed.

This script checks these two hypothesis via a sampling approach.
Perfect CL: Stage 1 labeling -> Mix nuclei -> Stage 2 labeling -> Fragment extraction / sequencing
Imperfect CL: Stage 1 labeling -> Fragment Extraction -> Mix fragments -> Stage 2 labeling -> sequencing
"""
import numpy as np
import pandas as pd

# Number of cells and wells in the first stage
cells_per_well_stage1 = 2500
labels_stage1 = 48
fragment_length = 250
genome_size = 175000000
# this might be much lower
fragments_per_cell = genome_size / fragment_length
fragments_per_cell = 500


def perfect_crosslinking_case():
    # Label all nuclei for all wells
    indices_stage1 = np.asarray([[str(i)] * cells_per_well_stage1
                                for i in range(labels_stage1)]).flatten()

    # Shuffle the nuclei
    np.random.shuffle(indices_stage1)

    indices_stage2 = []
    # take the indices from the previous round and append the new indices
    for second_index, start_pos in enumerate(
            range(0, len(indices_stage1), cells_per_well_stage1)):
        indices_stage2.append([indices_stage1[start_pos + x] +
                               '-' + str(second_index) for x in range(cells_per_well_stage1)])
    # finally, we extract the fragments_per_cell
    frags = [[ind] * fragments_per_cell for ind in indices_stage2]
    data = pd.DataFrame({'labels':np.asarray(frags).flatten()})
    return data


# Second, simulate the failed cross-linking scenario:

def no_crosslinking_case():
    # list of individual indices after the first stage
    # goes over cells and fragments
    indices_stage1 = np.asarray([[str(i)] * cells_per_well_stage1 * fragments_per_cell
                                  for i in range(labels_stage1)]).flatten()

    # First we simulate the under the perfect cross-linking scenario:
    np.random.shuffle(indices_stage1)
    indices_stage2 = []
    # take the indices from the previous round and append the new indices
    for second_index, start_pos in enumerate(
            range(0, len(indices_stage1), cells_per_well_stage1 * fragments_per_cell)):
        indices_stage2.append([split_indices_stage1[start_pos + x] +
                               '-' + str(second_index) for x in range(cells_per_well_stage1 * fragments_per_cell)])
    data = pd.DataFrame({'labels':np.asarray(indices_stage2).flatten()})
    return data


fig, axes = plt.subplots(nrows=2)

hyp1 = perfect_crosslinking_case()['labels'].value_counts()
hyp2 = no_crosslinking_case()['labels'].value_counts()

hyp1.plot(kind='bar', ax=axes[0], title='Perfect Cross-linking')
plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
hyp2.plot(kind='bar', ax=axes[1], title='No Cross-linking')

plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
plt.show()
