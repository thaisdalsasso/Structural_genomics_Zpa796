import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator

# Inputs
alignment_file = "./fasta/G.12_mature_mafft-linsi.fasta"
tm_score_file = "/Users/dalsasso/Desktop/Posdoc/CAU/structural_analysis/tm-align/Zpa796_tm-scores_all-vs-all_summary.txt"

##########################################
# Load alignment

basename = os.path.splitext(os.path.basename(alignment_file))[0]
alignment = AlignIO.read(alignment_file, "fasta")

sequence_ids = [record.id for record in alignment]
num_sequences = len(sequence_ids)


##########################################
# Sequence identity (symmetric matrix)

calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Sanity check
assert list(distance_matrix.names) == sequence_ids, \
    "Order mismatch between FASTA and DistanceCalculator output."

# Convert distance to symmetric matrix
distance_matrix_values = np.zeros((num_sequences, num_sequences), dtype=float)
for i in range(num_sequences):
    for j in range(i, num_sequences):
        dij = distance_matrix[i, j]
        distance_matrix_values[i, j] = dij
        distance_matrix_values[j, i] = dij

# Convert distance to percent identity: Identity (%) = (1 - Distance) * 100
identity_matrix = (1.0 - distance_matrix_values) * 100.0

# Save identity matrix 
identity_df = pd.DataFrame(identity_matrix, index=sequence_ids, columns=sequence_ids)
matrix_filename = f"{basename}_identity_matrix.csv"
identity_df.to_csv(matrix_filename)

# Generate heatmap
mpl.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(10, 8))
ax = sns.heatmap(
    identity_df, annot=False, cmap="Blues", cbar=True, vmin=0, vmax=100,
    cbar_kws={'label': 'Sequence Identity (%)'}
)
ax.set_title("")
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=18)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=18)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=22)
cbar.set_label('Sequence Identity (%)', fontsize=18)

heatmap_filename = f"{basename}_identity_heatmap.png"
plt.savefig(heatmap_filename, dpi=500, bbox_inches='tight')
plt.close()


##########################################
# TM-score (asymmetric matrix)

# Read TM-scores from file
tm_scores = {}
with open(tm_score_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 3:
            continue  
        query, subject, score = parts[:3]
        try:
            score_val = float(score)
        except ValueError:
            continue  
        tm_scores.setdefault(query, {})[subject] = score_val

# Warn if some alignment IDs never appear in the TM file at all
tm_ids = set(tm_scores.keys()) | {s for d in tm_scores.values() for s in d}
missing_from_tm = [sid for sid in sequence_ids if sid not in tm_ids]
if missing_from_tm:
    print("Warning: No TM-scores found for the following IDs (entire rows/cols may be NaN):")
    print("  " + ", ".join(missing_from_tm))

# Create asymmetric TM-score matrix
# Matrix[i, j] = TM-score(query_i -> subject_j) normalized by query_i length
tm_score_matrix = np.full((num_sequences, num_sequences), np.nan, dtype=float)

for i, query in enumerate(sequence_ids):
    for j, subject in enumerate(sequence_ids):
        if query == subject:
            tm_score_matrix[i, j] = 1.0  # Self-comparison
        else:
            tm_score_matrix[i, j] = tm_scores.get(query, {}).get(subject, np.nan)

# Report missing data (off-diagonal only)
off_diagonal_mask = ~np.eye(num_sequences, dtype=bool)
missing_count = int(np.isnan(tm_score_matrix[off_diagonal_mask]).sum())
total_comparisons = num_sequences * (num_sequences - 1)
missing_percentage = (missing_count / total_comparisons) * 100.0
print(f"Missing TM-score data: {missing_count}/{total_comparisons} ({missing_percentage:.1f}%)")

# Save TM-score matrix
tm_score_df = pd.DataFrame(tm_score_matrix, index=sequence_ids, columns=sequence_ids)
tm_score_filename = f"{basename}_tm_score_matrix_asymmetric.csv"
tm_score_df.to_csv(tm_score_filename)

# TM-score heatmap
plt.figure(figsize=(10, 8))
ax = sns.heatmap(
    tm_score_df, annot=False, cmap="Blues", cbar=True, vmin=0, vmax=1,
    mask=np.isnan(tm_score_matrix), cbar_kws={'label': 'TM-score'}
)
ax.set_title("")
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=18)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=18)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)
cbar.set_label('TM-score', fontsize=22)

tm_heatmap_filename = f"{basename}_tm_score_heatmap_asymmetric.png"
plt.savefig(tm_heatmap_filename, dpi=500, bbox_inches='tight')
plt.close()


##########################################
# Summary

print(f"Total sequences: {num_sequences}")
print(f"\nGenerated files:")
print(f"  1. {matrix_filename} (symmetric)")
print(f"  2. {heatmap_filename}")
print(f"  3. {tm_score_filename} (asymmetric)")
print(f"  4. {tm_heatmap_filename}")
print("\nNote: TM-score matrix is asymmetric because TM(A→B) ≠ TM(B→A).")
print("      Each score is normalized by the query protein's length.")
print("\nDone!")

