import os
import re

# Folder containing the TM-align results
folder_path = '/gxfs_work/cau/sunbo511/tmalign/Zpa796_all-vs-all'

# Output file for the summarized TM-scores
output_file = '/gxfs_work/cau/sunbo511/tmalign/Zpa796_secretome_tm-scores_all-vs-all_summary.txt'


# Extract protein names and TM-score from a file
def extract_tm_score(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        
        filename = os.path.basename(file_path)
        parts = filename.split('_vs_')
        if len(parts) != 2 or not parts[1].endswith('.tmalign.out'):
            raise ValueError(f"Filename {filename} does not match the expected pattern")
        protein1 = parts[0]
        protein2 = parts[1].replace('.tmalign.out', '')
        
        # Extract the first TM-score
        tm_score_match = re.search(r'TM-score=\s*(\d+\.\d+)', content)
        if not tm_score_match:
            raise ValueError(f"TM-score not found in file {file_path}")
        tm_score = float(tm_score_match.group(1))
        
        return protein1, protein2, tm_score


# Iterate over all files in the folder and extract TM-scores
with open(output_file, 'w') as out_file:
    out_file.write("#Protein1 (query)\tProtein2 (target)\tTM-score\n") # Write header
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.tmalign.out'):
                file_path = os.path.join(root, file)
                try:
                    protein1, protein2, tm_score = extract_tm_score(file_path)
                    out_file.write(f"{protein1}\t{protein2}\t{tm_score:.5f}\n")
                except Exception as e:
                    print(f"Error processing file {file_path}: {e}")



