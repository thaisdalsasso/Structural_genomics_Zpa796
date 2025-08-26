import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import PDB
from Bio.PDB.DSSP import DSSP
from Bio.SeqUtils import seq1

####################################################################
# Functions

# Calculate the Euclidean distance between two 3D coordinates
def calculate_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

# Extract C-alpha coordinates and residue names for a reference structure.
def get_ca_coordinates_and_residues(structure, chain_id):
    ca_coords = []
    residues = []
    for model in structure:
        chain = model[chain_id]
        model_ca_coords = []
        model_residues = []
        for residue in chain:
            if 'CA' in residue:
                model_ca_coords.append(residue['CA'].get_coord())
                residue_name = residue.get_resname()
                try:
                    one_letter_res = seq1(residue_name) 
                except KeyError:
                    one_letter_res = 'X'  # Unknown or non-standard residues
                model_residues.append(one_letter_res)
        ca_coords.append(np.array(model_ca_coords))
        residues.append(model_residues)
    return ca_coords, residues

# Extract protein names from PDB file
def get_model_ids_and_names(pdb_file_path):
    model_names = []
    with open(pdb_file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('REMARK') and 'Name:' in line:
                name = line.split('Name: ')[-1].strip()
                model_names.append(name)
    return model_names

# Calculate the secondary structure
def calculate_secondary_structure_dssp(structure, pdb_file_path):
    dssp = DSSP(structure[0], pdb_file_path)
    sec_structures = []
    
    for residue in dssp:
        ss = residue[2]
        if ss == 'H':
            sec_structures.append('a-helix')
        elif ss == 'E':
            sec_structures.append('b-strand')
        else:
            sec_structures.append('other')
    
    return sec_structures

# Calculate disulfide bonds
def get_disulfide_bonds(structure, chain_id):
    disulfide_bonds = []
    reference_chain = structure[0][chain_id]  # Use first model as reference

    cysteines = []
    for residue in reference_chain:
        res_name = residue.get_resname()
        res_id = residue.id[1]

        if res_name == "CYS" and residue.has_id("CA"):
            cysteines.append(residue)

    # Pair up cysteines and check distances based on Cα
    candidate_bonds = []
    for i, cys1 in enumerate(cysteines):
        for cys2 in cysteines[i+1:]:
            if abs(cys1.id[1] - cys2.id[1]) == 1:  # Prevent consecutive cysteines from bonding
                continue

            distance = calculate_distance(cys1["CA"].get_coord(), cys2["CA"].get_coord()) 
            if 2.0 <= distance <= 7.0:  # Valid disulfide bond range
                candidate_bonds.append((cys1.id[1], cys2.id[1], distance))

    candidate_bonds.sort(key=lambda x: x[2])
    assigned_cysteines = set()

    for cys1, cys2, dist in candidate_bonds:
        if cys1 in assigned_cysteines or cys2 in assigned_cysteines:
            continue  
        disulfide_bonds.append((cys1, cys2))
        assigned_cysteines.add(cys1)
        assigned_cysteines.add(cys2)

        print(f"Disulfide bond: Cys{cys1} - Cys{cys2}, {dist:.2f} Å")

    if not disulfide_bonds:
        print("No disulfide bonds detected.")

    return disulfide_bonds


# Generate line plot comparing each protein to a reference.
def plot_difference_to_reference(distance_matrix, residue_indices, sorted_protein_names, sec_structures, reference_protein_id, output_file, disulfide_bonds):
    plt.rcParams["font.family"] = "Arial"
    sns.set_style("white")

    color_palette = {
        0: "#0072B2",  # Blue
        1: "#D55E00",  # Vermillion
        2: "#F0E442",  # Yellow
        3: "#009E73",  # Green
        4: "#CC79A7",  # Pink
        5: "#56B4E9",  # Sky blue
        6: "#E69F00",  # Orange
        7: "#000000",  # Black
        8: "#999999",  # Gray
        9: "#AD7BE9",  # Lavender
        10: "#B79F00",  # Olive
        11: "#8172B2",  # Purple
        12: "#D5A676",  # Tan
        13: "#8C564B",  # Brown
        14: "#E377C2",  # Rose
    }

    # Calculate differences
    reference_index = sorted_protein_names.index(reference_protein_id)
    reference_distances = distance_matrix[reference_index, :]
    differences = distance_matrix - reference_distances

    plt.figure(figsize=(20, 8))

    # Secondary structure
    for i, ss in enumerate(sec_structures):
        if ss == 'a-helix':
            plt.axvspan(i, i + 1, color='lightblue', alpha=0.2)
        elif ss == 'b-strand':
            plt.axvspan(i, i + 1, color='lightsalmon', alpha=0.2)

    # Plot distances for all proteins
    for i, protein_name in enumerate(sorted_protein_names):
        if protein_name != reference_protein_id: 
            plt.plot(
                residue_indices, 
                differences[i, :], 
                label=f'{protein_name}', 
                color=color_palette[i % 15],
                alpha=0.8,  
                linewidth=3.0
            )

    #y_max = np.max(differences) + 10 
    y_max = 75
    y_offset_step = 5  
    current_y = y_max + 5  

    # Annotate disulfide bonds 
    for cys1, cys2 in disulfide_bonds:
        cys1_x = cys1 - 0.5  
        cys2_x = cys2 - 0.5  
        
        plt.plot([cys1_x, cys2_x], [current_y, current_y], color='black', linestyle='-', linewidth=1)  
        plt.scatter([cys1_x, cys2_x], [current_y, current_y], color='black', s=16, zorder=2)  
        plt.text((cys1_x + cys2_x) / 2, current_y + 1, f"Cys{cys1}-Cys{cys2}", 
                 ha='center', va='bottom', fontsize=14, color='black')

        current_y += y_offset_step 

    plt.xlabel("Residue index", fontsize=20)
    plt.ylabel("Distance (Å)", fontsize=20)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    x_ticks = [1] + list(np.arange(5, max(residue_indices) + 1, 5))
    plt.xticks([x - 0.5 for x in x_ticks], labels=[str(x) for x in x_ticks], fontsize=18)
    plt.yticks(np.arange(0, current_y + 10, 10), fontsize=18)
    plt.ylim(0, current_y + 5)  
    plt.legend(fontsize=20, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
    plt.savefig(output_file.replace('.csv', '_distance_differences_lineplot.svg'), dpi=500, bbox_inches='tight')
    plt.close()

    return color_palette


####################################################################
# Load PDB file and specify reference protein ID

pdb_file_path = './G.12_foldmason.pdb' 
reference_protein_id = "jg1.t1"   # The reference is assumed to be the shortest and first protein (model 0) in the alignment

####################################################################
# Apply functions and write outputs

pdb_parser = PDB.PDBParser()
structure = pdb_parser.get_structure("aligned", pdb_file_path)

ca_coords_list, residues_list = get_ca_coordinates_and_residues(structure, 'A')

protein_names = get_model_ids_and_names(pdb_file_path)

sorted_protein_names = sorted(protein_names)

sorted_ca_coords_list = [ca_coords_list[protein_names.index(name)] for name in sorted_protein_names]

# Find the index of the reference protein
try:
    reference_model_index = sorted_protein_names.index(reference_protein_id)
except ValueError:
    raise ValueError(f"Protein ID '{reference_protein_id}' not found in the PDB file.")

reference_coords = sorted_ca_coords_list[reference_model_index]
reference_residues = residues_list[reference_model_index]

disulfide_bonds = get_disulfide_bonds(structure, 'A')

basename = os.path.splitext(os.path.basename(pdb_file_path))[0]
output_file = f"{basename}_Calpha_euclidean_distance_secondary-structures.csv"

sec_structures = calculate_secondary_structure_dssp(structure, pdb_file_path)

distance_matrix = []
residue_indices = []

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    header = ["Residue Index", "Amino Acid", "Secondary Structure"] + sorted_protein_names
    writer.writerow(header)
    
    for i in range(len(reference_coords)):
        residue_distances = [i + 1, reference_residues[i], sec_structures[i]]
        current_distances = []  
        
        for model_coords in sorted_ca_coords_list:
            aligned_coord = model_coords[i]
            ref_coord = reference_coords[i]
            distance = calculate_distance(ref_coord, aligned_coord)
            residue_distances.append(distance)
            current_distances.append(distance)
        
        writer.writerow(residue_distances)
        distance_matrix.append(current_distances)
        residue_indices.append(i + 1)

distance_matrix = np.array(distance_matrix).T  

plot_difference_to_reference(distance_matrix, residue_indices, sorted_protein_names, sec_structures, reference_protein_id, output_file, disulfide_bonds)

