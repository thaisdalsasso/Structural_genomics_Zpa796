import pandas as pd
import requests
from bs4 import BeautifulSoup

def extract_scope_classification(scop_id):
    url = f'http://scop.berkeley.edu/sunid={scop_id}'
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')
    
    classification = {
        "Class": None,
        "Fold": None,
        "Superfamily": None,
        "Family": None
    }
    
    lineage = soup.find('div', {'class': 'col-md-12'}).find('ol', {'class': 'browse'})
    
    if lineage:
        items = lineage.find_all('li')
        for item in items:
            text = item.get_text(strip=True)
            if 'Class' in text:
                classification['Class'] = text.split('Class')[-1].split('[')[0].strip()
            elif 'Fold' in text:
                classification['Fold'] = text.split('Fold')[-1].split('[')[0].strip()
            elif 'Superfamily' in text:
                classification['Superfamily'] = text.split('Superfamily')[-1].split('[')[0].strip()
            elif 'Family' in text:
                classification['Family'] = text.split('Family')[-1].split('[')[0].strip()
    
    return classification

# Function to retrieve the SCOPe family ID from the table
def retrieve_scope_id(filename, target_id):
    with open(filename, 'r') as file:
        for line in file:
            columns = line.split()
            if columns[0] == target_id:
                fa_id = [part.split('=')[1] for part in columns[-1].split(',') if part.startswith('fa=')]
                if fa_id:
                    return fa_id[0]
    return None

# Load the best hits data
best_hits_path = './easy-search/Zpa796_AF2bestmodels.vs.SCOPe40_foldseek_qcov-qtmscore-besthits.tsv'
best_hits_data = pd.read_csv(best_hits_path, sep='\t')
print("Loaded best hits data.")

# Extract SCOPe data for each entry in the best hits
def extract_scope_annotations(row, filename):
    target_id = row['target'].split('.ent')[0]
    scope_id = retrieve_scope_id(filename, target_id)
    print(f"Retrieving data for: {scope_id}")
    
    if scope_id:
        scope_info = extract_scope_classification(scope_id)
        return pd.Series({
            'scope_class': scope_info.get('Class', ''),
            'scope_fold': scope_info.get('Fold', ''),
            'scope_superfamily': scope_info.get('Superfamily', ''),
            'scope_family': scope_info.get('Family', '')
        })
    else:
        return pd.Series({
            'scope_class': '',
            'scope_fold': '',
            'scope_superfamily': '',
            'scope_family': ''
        })

# Filepath for SCOPe table
scope_table_path = '/Users/dalsasso/foldseek/info_databases/SCOPe_v2.08/SCOPe_F40_2.08_descriptions/dir.cla.scope.2.08-stable.txt'

annotations = best_hits_data.apply(lambda row: extract_scope_annotations(row, scope_table_path), axis=1)
best_hits_data = pd.concat([best_hits_data, annotations], axis=1)
print("SCOPe annotations extracted.")

# Save to a new TSV file
output_path = './easy-search/Zpa796_AF2bestmodels.vs.SCOPe40_foldseek_qcov-qtmscore-besthits_entry-anno.tsv'
best_hits_data.to_csv(output_path, sep='\t', index=False)
print(f"Annotated data saved to {output_path}.")

