# Rebuilds the 61-channel adjacency used by the evaluation scripts.
# Step 1: MNE elec1010 template, restricted to the 61 channels (190 pairs).
# Step 2: ten extra pairs for sparsely connected edge electrodes (200 pairs in total).
# Writes adjacency_61ch.csv; adjacency_csv_to_rds.R converts it to adjacency_xproj_yproj_61ch.rds.
import mne
import numpy as np
import pandas as pd

ch = ['AF1', 'AF2', 'AF7', 'AF8', 'AFZ', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6',
      'CP1', 'CP2', 'CP3', 'CP4', 'CP5', 'CP6', 'CPZ', 'CZ',
      'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8',
      'FC1', 'FC2', 'FC3', 'FC4', 'FC5', 'FC6', 'FCZ',
      'FP1', 'FP2', 'FPZ', 'FT7', 'FT8', 'FZ',
      'O1', 'O2', 'OZ', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8',
      'PO1', 'PO2', 'PO7', 'PO8', 'POZ', 'PZ',
      'T7', 'T8', 'TP7', 'TP8']
extra = [('AF1', 'AF7'), ('AF2', 'AF8'), ('AF7', 'F3'), ('AF8', 'F4'), ('P3', 'PO7'),
         ('PO1', 'PO7'), ('P4', 'PO8'), ('PO2', 'PO8'), ('O1', 'POZ'), ('O2', 'POZ')]

adj, names = mne.channels.read_ch_adjacency('elec1010')
dense = adj.toarray().astype(int)
pos = {n.upper(): i for i, n in enumerate(names)}
idx = [pos[c] for c in ch]
sub = dense[np.ix_(idx, idx)]
np.fill_diagonal(sub, 0)
assert sub.sum() // 2 == 190
df = pd.DataFrame(sub, index=ch, columns=ch)
for a, b in extra:
    df.loc[a, b] = 1
    df.loc[b, a] = 1
assert df.values.sum() // 2 == 200
df.to_csv('adjacency_61ch.csv')
print('MNE', mne.__version__, 'pairs', df.values.sum() // 2)
