
import pandas as pd
import os


input_file_path = "./easy-search/Zpa796_AF2bestmodels.vs.SCOPe40_foldseek.tsv"


##################################################################
# Prepare data

data = pd.read_csv(input_file_path, sep='\t')
data['qcov'] = pd.to_numeric(data['qcov'], errors='coerce')
data['qtmscore'] = pd.to_numeric(data['qtmscore'], errors='coerce')

# Filter by max qcov
max_qcov_per_query = data.groupby('query')['qcov'].max().reset_index()
max_qcov_per_query.rename(columns={'qcov': 'max_qcov'}, inplace=True)

data_merged = pd.merge(data, max_qcov_per_query, on='query')
best_qcov_hits = data_merged[data_merged['qcov'] == data_merged['max_qcov']]

# From those with the highest qcov, select the one with the highest qtmscore
best_hits = best_qcov_hits.loc[best_qcov_hits.groupby('query')['qtmscore'].idxmax()]


# Output file 
base_name = os.path.splitext(os.path.basename(input_file_path))[0]  
output_file_name = f"{base_name}_qcov-qtmscore-besthits.tsv"
output_file_path = os.path.join(os.path.dirname(input_file_path), output_file_name)

best_hits.to_csv(output_file_path, sep='\t', index=False)

print(f"Number of unique best entries: {len(best_hits)}")
