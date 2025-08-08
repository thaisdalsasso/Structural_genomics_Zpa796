import os
import shutil


# Directory with AF2 results
source_folder = "./Zpa796_AF2"

# Directpry for AF2 best models only
destination_folder = "./Zpa796_AF2bestmodels"


os.makedirs(destination_folder, exist_ok=True)


for protein_id in os.listdir(source_folder):
    protein_path = os.path.join(source_folder, protein_id)
    
    if not os.path.isdir(protein_path):
        continue  # Skip files

    ranked_pdb_file = os.path.join(protein_path, "ranked_0.pdb")
    if os.path.exists(ranked_pdb_file):
        destination_file = os.path.join(destination_folder, f"{protein_id}.pdb")
        shutil.copy(ranked_pdb_file, destination_file)
        print(f"Copied: {protein_id}/ranked_0.pdb -> {destination_file}")
    else:
        print(f"ranked_0.pdb not found in {protein_path}")

print("Done.")

