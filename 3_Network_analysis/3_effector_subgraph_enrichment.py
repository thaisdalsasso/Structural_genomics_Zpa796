import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt


# File names:

file_path = '/Users/dalsasso/Desktop/Posdoc/CAU/annotations/Zpa796_secretome_metadata.tsv'  
output_file_path = './effector_subgraph_enrichment_results.tsv'
output_plot_file = './effector_subgraph_enrichment_cluster_size.png'


##############################################################################################################
# Enrichment test:

data = pd.read_csv(file_path, sep='\t')

def map_effector_status(category):
    if pd.isna(category):
        return category  # Leave NaN as is
    elif 'Non-effector' in category:
        return 'Non-effector'
    else:
        return 'Effector'

data['EffectorP3'] = data['EffectorP3'].apply(map_effector_status)

filtered_data = data[['Structural subgraph', 'EffectorP3']].dropna()

# Create a contingency table
contingency_table = pd.crosstab(filtered_data['Structural subgraph'], filtered_data['EffectorP3'])

# Calculate Fisher's exact test
results = []
p_values = []
cluster_sizes = []
for index, row in contingency_table.iterrows():
    if len(row) == 2:  # Ensuring we have a complete row for the test
        table = [row.tolist(), [contingency_table.iloc[:,0].sum() - row.iloc[0], contingency_table.iloc[:,1].sum() - row.iloc[1]]]
        #odds_ratio, p_value = fisher_exact(table, alternative='two-sided') # test for both under/overrepresented
        odds_ratio, p_value = fisher_exact(table, alternative='greater') #test only for overrepresented
        results.append((index, odds_ratio, p_value))
        p_values.append(p_value)
        cluster_sizes.append(sum(row))  

'''
Contingency table structure per cluster:
-------------------------------------------------------
                          | Effector | Non-effector
-------------------------------------------------------
In this cluster           |    a     |       b
Outside this cluster      |    c     |       d

Where:
- a = Number of effectors in the cluster (row.iloc[0])
- b = Number of non-effectors in the cluster (row.iloc[1])
- c = Total number of effectors outside the cluster = total_effectors - a
- d = Total number of non-effectors outside the cluster = total_non_effectors - b
'''


# Apply correction
corrected_p_values = multipletests(p_values, alpha=0.05, method='fdr_bh')[1]

corrected_results = [(*res, cpv, size) for res, cpv, size in zip(results, corrected_p_values, cluster_sizes)]

results_df = pd.DataFrame(corrected_results, columns=['Structural subgraph', 'Odds Ratio', 'P-value', 'Adjusted P-value', 'Cluster Size'])
results_df.to_csv(output_file_path, sep='\t', index=False)

significant_results = results_df[results_df['Adjusted P-value'] < 0.05].sort_values(by='Adjusted P-value')
print(significant_results)


##############################################################################################################
# Generate cluster size plot

plt.rcParams["font.family"] = "Arial"

plt.figure(figsize=(10, 5))

# Raw p-values
plt.subplot(1, 2, 1)
plt.scatter(results_df['Cluster Size'], results_df['P-value'], color='steelblue', alpha=0.5)  
plt.axhline(y=0.05, color='r', linestyle='--')
plt.title('Cluster Size vs. Raw P-value')
plt.xlabel('Cluster Size')
plt.ylabel('P-value')

# Adjusted p-values
plt.subplot(1, 2, 2)
plt.scatter(results_df['Cluster Size'], results_df['Adjusted P-value'], color='teal', alpha=0.5) 
plt.axhline(y=0.05, color='r', linestyle='--')
plt.title('Cluster Size vs. Adjusted P-value')
plt.xlabel('Cluster Size')
plt.ylabel('Adjusted P-value')

plt.tight_layout()

# Save plot
plt.savefig(output_plot_file, format='png', dpi=300) 

print(f"Plot saved as {output_plot_file}")

