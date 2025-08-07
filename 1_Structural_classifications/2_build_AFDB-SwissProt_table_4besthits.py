import pandas as pd
import requests

def extract_uniprot_id(filename):
    return filename.split('-')[1]

def get_uniprot_data(uniprot_id):
    print(f"Retrieving UniProt data for ID: {uniprot_id}")
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        entry_name = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'No entry name available')
        molecular_functions = ', '.join([keyword['name'] for keyword in data.get('keywords', []) if keyword.get('category') == 'Molecular function'])
        return {
            'uniprot_entry_name': entry_name,
            'uniprot_molecular_function': molecular_functions
        }
    else:
        print(f"Failed to retrieve UniProt data for {uniprot_id}, status code: {response.status_code}")
        return {'uniprot_entry_name': 'No entry name available', 'uniprot_molecular_function': ''}

# Load the best hits data
best_hits_path = './easy-search/Zpa796_AF2bestmodels.vs.AFDB-SwissProt_foldseek_qcov-qtmscore-besthits.tsv'
best_hits_data = pd.read_csv(best_hits_path, sep='\t')
print("Loaded best hits data.")

# Extract UniProt annotations for each entry in the best hits
def extract_uniprot_annotations(row):
    uniprot_id = extract_uniprot_id(row['target'])
    uniprot_info = get_uniprot_data(uniprot_id)
    return pd.Series({
        'uniprot_id': uniprot_id,
        'uniprot_entry_name': uniprot_info['uniprot_entry_name'],
        'uniprot_molecular_function': uniprot_info['uniprot_molecular_function']
    })

uniprot_annotations = best_hits_data.apply(extract_uniprot_annotations, axis=1)
best_hits_data = pd.concat([best_hits_data, uniprot_annotations], axis=1)
print("UniProt annotations extracted.")

# Save to a new TSV file
output_path = './easy-search/Zpa796_AF2bestmodels.vs.AFDB-SwissProt_foldseek_qcov-qtmscore-besthits_uniprot-anno.tsv'
best_hits_data.to_csv(output_path, sep='\t', index=False)
print(f"Annotated data saved to {output_path}.")

