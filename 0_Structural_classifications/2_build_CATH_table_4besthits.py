import pandas as pd
import requests
from bs4 import BeautifulSoup


def get_cath_data_by_id(cath_id):
    if cath_id.startswith('af_'):
        cath_code = cath_id.split('_')[-1]
        return get_cath_data_from_code(cath_code)
    else:
        return get_cath_data_from_cath_id(cath_id)

def get_cath_data_from_cath_id(cath_id):
    url = f"https://www.cathdb.info/version/latest/domain/{cath_id}"
    print(f"Retrieving data from CATH domain for ID: {cath_id}")
    response = requests.get(url)
    if response.status_code == 200:
        return parse_cath_response(response.content)
    else:
        print(f"Failed to retrieve data for {cath_id}, status code: {response.status_code}")
        return {}

def get_cath_data_from_code(cath_code):
    url = f"https://www.cathdb.info/version/latest/superfamily/{cath_code}/classification"
    print(f"Retrieving classification from CATH for code: {cath_code}")
    response = requests.get(url)
    if response.status_code == 200:
        return parse_cath_response(response.content)
    else:
        print(f"Failed to retrieve classification for {cath_code}, status code: {response.status_code}")
        return {}

def parse_cath_response(content):
    print("Parsing CATH response.")
    soup = BeautifulSoup(content, 'html.parser')
    table = soup.find('table', class_='table-condensed')
    data = {}
    if table:
        rows = table.find_all('tr')[1:]  # Skip header row
        for row in rows:
            columns = row.find_all('td')
            if len(columns) >= 3:
                level = columns[0].img['alt'] if columns[0].find('img') else None
                code = columns[1].get_text(strip=True)
                description = columns[2].get_text(strip=True)
                if level:
                    data[level] = f"{description} ({code})"
    return data

# Load the best hits data
best_hits_path = './easy-search/Zpa796_AF2bestmodels.vs.CATH50_foldseek_qcov-qtmscore-besthits.tsv'
best_hits_data = pd.read_csv(best_hits_path, sep='\t')
print("Loaded best hits data.")

# Extract CATH data for each entry in the best hits
def extract_cath_annotations(row):
    print(f"Processing {row['target']}")
    cath_info = get_cath_data_by_id(row['target'])
    return pd.Series({
        'cath_class (code)': cath_info.get('Class', ''),
        'cath_architecture (code)': cath_info.get('Architecture', ''),
        'cath_topology (code)': cath_info.get('Topology', ''),
        'cath_superfamily (code)': cath_info.get('Homologous Superfamily', '')
    })

annotations = best_hits_data.apply(extract_cath_annotations, axis=1)
best_hits_data = pd.concat([best_hits_data, annotations], axis=1)
print("CATH annotations extracted.")

# Save to a new TSV file
output_path = './easy-search/Zpa796_AF2bestmodels.vs.CATH50_foldseek_qcov-qtmscore-besthits_entry-anno.tsv'
best_hits_data.to_csv(output_path, sep='\t', index=False)
print(f"Annotated data saved to {output_path}.")
