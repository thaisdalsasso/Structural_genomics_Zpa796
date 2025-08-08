import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

mpl.rcParams['font.family'] = 'Arial'

# Define custom colors for annotations
hydro_cmap = LinearSegmentedColormap.from_list("hydro", ["darkcyan", "white", "darkgoldenrod"])

# Hydrophobicity 
kyte_doolittle_scale = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}


def parse_coordinates_and_sequences(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    aligned_atoms = []
    sequence = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    aligned_atoms.append(residue['CA'].get_coord())
                    aa = residue.resname.capitalize()[0]  
                    sequence.append(aa)

    return np.array(aligned_atoms), sequence

def compute_distance_matrix(coords1, coords2):
    return np.sqrt(np.sum((coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]) ** 2, axis=2))

def compute_hydrophobicity(sequence):
    hydrophobicity = [kyte_doolittle_scale.get(aa, 0) for aa in sequence]
    return np.array(hydrophobicity)

def extract_protein_ids(filename):
    base_name = os.path.basename(filename).replace(".pdb", "")
    protein1, protein2 = base_name.split("_vs_")
    return protein1, protein2, base_name

def plot_contact_map_with_annotations(distances, protein1, protein2, 
                                      protein1_hydro, protein2_hydro, base_name):

    fig, ax = plt.subplots(figsize=(12, 12))

    # Debugging info
    print(f"Protein 1: {protein1}, Length: {len(protein1_hydro)}")
    print(f"Protein 2: {protein2}, Length: {len(protein2_hydro)}")
    print(f"Hydrophobicity Protein 1: {protein1_hydro}")
    print(f"Hydrophobicity Protein 2: {protein2_hydro}")

    # Heatmap for distance matrix
    sns.heatmap(
        distances, cmap="rocket", cbar=False, vmax=10, ax=ax
    )

    hydro_vmin, hydro_vmax = -4, 4

    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_xlabel(f"{protein2}", fontsize=32, fontweight="bold", labelpad=10)
    ax.set_ylabel(f"{protein1}", fontsize=32, fontweight="bold", labelpad=10)

    ax.text(-1, -6, "N", fontsize=28, ha="right", va="bottom", color="black")  # N on y-axis
    #ax.text(-4, -4, "N", fontsize=28, ha="center", va="bottom", color="black")  # N on x-axis
    ax.text(-1 , len(protein1_hydro) -3, "C", fontsize=28, ha="right", va="center", color="black")  # C on y-axis
    ax.text(len(protein2_hydro) - 3, -6, "C", fontsize=28, ha="left", va="bottom", color="black")  # C on x-axis

    # Hydrophobicity annotation
    ax_hydro_x = fig.add_axes([ax.get_position().x0, ax.get_position().y1 + 0.01, ax.get_position().width, 0.05])
    ax_hydro_y = fig.add_axes([ax.get_position().x1 + 0.01, ax.get_position().y0, 0.05, ax.get_position().height])

    ax_hydro_x.imshow(protein2_hydro.reshape(1, -1), aspect="auto", cmap=hydro_cmap, vmin=hydro_vmin, vmax=hydro_vmax)
    ax_hydro_y.imshow(protein1_hydro[::-1].reshape(-1, 1), aspect="auto", cmap=hydro_cmap, vmin=hydro_vmin, vmax=hydro_vmax)

    ax_hydro_x.set_xticks([]); ax_hydro_x.set_yticks([])
    ax_hydro_y.set_xticks([]); ax_hydro_y.set_yticks([])


    # Legends
    legend_fontsize = 34  
    ticks_fontsize = 28   
    
    cbar_ax_hydro = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.12, ax.get_position().width, 0.03])
    cbar_ax_dist = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.24, ax.get_position().width, 0.03])
    
    cbar_hydro = fig.colorbar(plt.cm.ScalarMappable(cmap=hydro_cmap, norm=plt.Normalize(vmin=hydro_vmin, vmax=hydro_vmax)), cax=cbar_ax_hydro, orientation='horizontal')
    cbar_dist = fig.colorbar(plt.cm.ScalarMappable(cmap="rocket", norm=plt.Normalize(vmin=0, vmax=10)), cax=cbar_ax_dist, orientation='horizontal')
    
    cbar_hydro.set_label("Hydrophobicity", fontsize=legend_fontsize)
    cbar_dist.set_label("Distance (Å)", fontsize=legend_fontsize)

    cbar_hydro.ax.tick_params(labelsize=ticks_fontsize)
    cbar_dist.ax.tick_params(labelsize=ticks_fontsize)

    plt.savefig(f"{base_name}_contact_map_annotated.png", dpi=300, bbox_inches="tight")
    #plt.show()

# Input alignment and execute functions
pdb_file = "jg9741.t1_vs_ZtIPO323_014270.1.pdb" # Pairwise alignment from TM-align
protein1, protein2, base_name = extract_protein_ids(pdb_file)

coordinates, sequences = parse_coordinates_and_sequences(pdb_file)

split_index = len(coordinates) // 2
coords1, coords2 = coordinates[:split_index], coordinates[split_index:]

hydro1 = compute_hydrophobicity(sequences[:split_index])
hydro2 = compute_hydrophobicity(sequences[split_index:])

plot_contact_map_with_annotations(compute_distance_matrix(coords1, coords2), protein1, protein2, hydro1, hydro2, base_name)

