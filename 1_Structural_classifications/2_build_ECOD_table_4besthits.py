import pandas as pd
import requests
from bs4 import BeautifulSoup

def extract_ecod_classification(ecod_id):
    url = f"http://prodata.swmed.edu/ecod/complete/domain/{ecod_id}"
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Failed to retrieve data for ECOD ID: {ecod_id}")
        return None
    
    soup = BeautifulSoup(response.content, 'html.parser')
    
    classification = {
        "A-group": None,
        "X-group": None,
        "H-group": None,
        "T-group": None,
        "F-group": None
    }
    
    try:
        li_elements = soup.find('div', class_='main-panel').find_all('li')
        for li in li_elements:
            text = li.get_text(strip=True)
            if text.startswith('A:'):
                classification['A-group'] = text.split(': ')[-1]
            elif text.startswith('X:'):
                classification['X-group'] = text.split(': ')[-1]
            elif text.startswith('H:'):
                classification['H-group'] = text.split(': ')[-1]
            elif text.startswith('T:'):
                classification['T-group'] = text.split(': ')[-1]
            elif text.startswith('F:'):
                classification['F-group'] = text.split(': ')[-1]
    except AttributeError as e:
        print(f"Error parsing ECOD classification: {e}")
    
    return classification

# Function to retrieve the ECOD ID from the table
def retrieve_ecod_id_from_table(file_path, search_id):
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.split('\t')  
            if columns[0] == search_id:
                ecod_id = columns[1]  
                #print(f"Retrieved ECOD ID: {ecod_id} for Search ID: {search_id}")
                return ecod_id
    return None



# Load the best hits data
best_hits_path = './easy-search/Zpa796_AF2bestmodels.vs.ECOD40_foldseek_qcov-qtmscore-besthits.tsv'
best_hits_data = pd.read_csv(best_hits_path, sep='\t')
print("Loaded best hits data.")

# Extract ECOD data for each entry in the best hits
def extract_ecod_annotations(row, filename):
    search_id = row['target'].split('.pdbnum')[0]  
    ecod_id = retrieve_ecod_id_from_table(filename, search_id)
    print(f"Retrieving data for: {ecod_id}")

    if ecod_id:
        ecod_info = extract_ecod_classification(ecod_id)
        return pd.Series({
            'ecod_id': ecod_id,
            'ecod_a_group': ecod_info.get('A-group', ''),
            'ecod_x_group': ecod_info.get('X-group', ''),
            'ecod_h_group': ecod_info.get('H-group', ''),
            'ecod_t_group': ecod_info.get('T-group', ''),
            'ecod_f_group': ecod_info.get('F-group', '')
        })
    else:
        return pd.Series({
            'ecod_id': '',
            'ecod_a_group': '',
            'ecod_x_group': '',
            'ecod_h_group': '',
            'ecod_t_group': '',
            'ecod_f_group': ''
        })

# Filepath for ECOD table
ecod_table_path = '/Users/dalsasso/foldseek/info_databases/ECOD_v291/ecod.latest.domains.txt'

annotations = best_hits_data.apply(lambda row: extract_ecod_annotations(row, ecod_table_path), axis=1)
best_hits_data = pd.concat([best_hits_data, annotations], axis=1)
print("ECOD annotations extracted.")

# Save to a new TSV file
output_path = './easy-search/Zpa796_AF2bestmodels.vs.ECOD40_foldseek_qcov-qtmscore-besthits_entry-anno.tsv'
best_hits_data.to_csv(output_path, sep='\t', index=False)
print(f"Annotated data saved to {output_path}.")
