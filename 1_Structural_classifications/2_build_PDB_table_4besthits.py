import pandas as pd
import requests

def extract_pdb_id(filename):
    return filename.split('-')[0]

def get_pdb_data(pdb_id):
    print(f"Retrieving PDB data for ID: {pdb_id}")
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        description = data.get('struct', {}).get('title', 'No description available')
        classifications = data.get('struct_keywords', {}).get('pdbx_keywords', 'No classification available')
        return {
            'pdb_description': description,
            'pdb_classification': classifications
        }
    else:
        print(f"Failed to retrieve PDB data for {pdb_id}, status code: {response.status_code}")
        return {'pdb_description': 'No description available', 'pdb_classification': 'No classification available'}

# Load the best hits data
best_hits_path = './easy-search/Zpa796_AF2bestmodels.vs.PDB_foldseek_qcov-qtmscore-besthits.tsv'
best_hits_data = pd.read_csv(best_hits_path, sep='\t')
print("Loaded best hits data.")

# Extract PDB annotations for each entry in the best hits
def extract_pdb_annotations(row):
    pdb_id = extract_pdb_id(row['target'])
    pdb_info = get_pdb_data(pdb_id)
    return pd.Series({
        'pdb_id': pdb_id,
        'pdb_description': pdb_info['pdb_description'],
        'pdb_classification': pdb_info['pdb_classification']
    })

pdb_annotations = best_hits_data.apply(extract_pdb_annotations, axis=1)
best_hits_data = pd.concat([best_hits_data, pdb_annotations], axis=1)
print("PDB annotations extracted.")

# Save to a new TSV file
output_path = './easy-search/Zpa796_AF2bestmodels.vs.PDB_foldseek_qcov-qtmscore-besthits_pdb-anno.tsv'
best_hits_data.to_csv(output_path, sep='\t', index=False)
print(f"Annotated data saved to {output_path}.")

