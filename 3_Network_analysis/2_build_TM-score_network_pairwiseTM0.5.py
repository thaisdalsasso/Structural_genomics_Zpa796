import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from kneed import KneeLocator
import numpy as np

# Define files:
file_path = 'Zpa796_secretome_tm-scores_all-vs-all_summary.txt'
output_file = 'Zpa796_secretome_tm-scores_network.graphml'
stats_output_file = 'Zpa796_secretome_tm-scores_network_stats.txt'
elbow_plot_file = 'Zpa796_secretome_tm-scores_network_elbow_plot.pdf'

# Function to read and filter data
def read_and_filter_data(file_path):
    data = {}
    with open(file_path, 'r') as file:
        next(file)
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                protein1 = parts[0]
                protein2 = parts[1]
                try:
                    tm_score = float(parts[2])
                    pair = tuple(sorted((protein1, protein2)))
                    if pair not in data:
                        data[pair] = []
                    data[pair].append(tm_score)
                except ValueError:
                    print(f"Invalid TM-score value: {parts[2]}")
            else:
                print(f"Invalid line format: {line}")
    return data


data = read_and_filter_data(file_path)

highest_scores = {}

# Filter bi-direction comparisons TM-score = 0.5 and get the highest value for each pair 
for pair, scores in data.items():
    if len(scores) == 2 and all(score >= 0.5 for score in scores):
        highest_scores[pair] = max(scores)

# Create a graph and add edges with the highest TM-score for each pair
G = nx.Graph()

for (protein1, protein2), tm_score in highest_scores.items():
    G.add_edge(protein1, protein2, weight=tm_score)

# Compute network statistics
degree_dict = dict(G.degree())
degrees = list(degree_dict.values())
average_degree = sum(degrees) / float(G.number_of_nodes())
density = nx.density(G)
average_clustering_coefficient = nx.average_clustering(G)

# Handle connected and disconnected graphs
if nx.is_connected(G):
    average_shortest_path_length = nx.average_shortest_path_length(G)
    diameter = nx.diameter(G)
else:
    # Find the largest connected component
    largest_cc = max(nx.connected_components(G), key=len)
    subgraph = G.subgraph(largest_cc)
    average_shortest_path_length = nx.average_shortest_path_length(subgraph)
    diameter = nx.diameter(subgraph)

connected_components = nx.number_connected_components(G)

# Compute centrality measures
degree_centrality = nx.degree_centrality(G)
betweenness_centrality = nx.betweenness_centrality(G)
closeness_centrality = nx.closeness_centrality(G)

avg_degree_centrality = sum(degree_centrality.values()) / len(degree_centrality)
avg_betweenness_centrality = sum(betweenness_centrality.values()) / len(betweenness_centrality)
avg_closeness_centrality = sum(closeness_centrality.values()) / len(closeness_centrality)


# Save network statistics
with open(stats_output_file, 'w') as f:
    f.write(f"Number of nodes: {G.number_of_nodes()}\n")
    f.write(f"Number of edges: {G.number_of_edges()}\n")
    f.write(f"Average degree: {average_degree:.3f}\n")
    f.write(f"Density: {density:.3f}\n")
    f.write(f"Average clustering coefficient: {average_clustering_coefficient:.3f}\n")
    f.write(f"Average shortest path length: {average_shortest_path_length:.3f}\n")
    f.write(f"Diameter: {diameter}\n")
    f.write(f"Number of connected components: {connected_components}\n")
    f.write(f"\nAverage degree centrality: {avg_degree_centrality:.3f}\n")
    f.write(f"Average betweenness centrality: {avg_betweenness_centrality:.3f}\n")
    f.write(f"Average closeness centrality: {avg_closeness_centrality:.3f}\n")

print(f"Network statistics saved to {stats_output_file}")

# Add centrality measures as node attributes
nx.set_node_attributes(G, {k: round(v, 3) for k, v in degree_centrality.items()}, 'degree_centrality')
nx.set_node_attributes(G, {k: round(v, 3) for k, v in betweenness_centrality.items()}, 'betweenness_centrality')
nx.set_node_attributes(G, {k: round(v, 3) for k, v in closeness_centrality.items()}, 'closeness_centrality')
# nx.set_node_attributes(G, {k: round(v, 3) for k, v in eigenvector_centrality.items()}, 'eigenvector_centrality')


# Perform K-means clustering and add cluster labels
k_range = range(1, 21)
def perform_kmeans_clustering(centrality_values, G, attribute_name):
    scaler = StandardScaler()
    centrality_values_scaled = scaler.fit_transform(np.array(centrality_values).reshape(-1, 1))

    # Elbow method to determine the optimal number of clusters
    sse = []
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(centrality_values_scaled)
        sse.append(kmeans.inertia_)

    kneedle = KneeLocator(k_range, sse, curve='convex', direction='decreasing')
    optimal_k = kneedle.elbow

    # K-means clustering
    kmeans = KMeans(n_clusters=optimal_k, random_state=42)
    kmeans.fit(centrality_values_scaled)

    # Compute the mean centrality values for each cluster
    centrality_means = []
    for cluster_label in range(optimal_k):
        cluster_indices = np.where(kmeans.labels_ == cluster_label)[0]
        mean_centrality = np.mean([centrality_values[i] for i in cluster_indices])
        centrality_means.append((cluster_label, mean_centrality))

    sorted_clusters = sorted(centrality_means, key=lambda x: x[1])

    cluster_mapping = {old_label: chr(65 + i) for i, (old_label, _) in enumerate(sorted_clusters)}

    # Add sorted cluster labels as node attributes
    for i, node in enumerate(G.nodes):
        old_cluster_label = kmeans.labels_[i]
        G.nodes[node][f'{attribute_name}_cluster'] = cluster_mapping[old_cluster_label]

    return kmeans, optimal_k, sse

degree_kmeans, degree_optimal_k, degree_sse = perform_kmeans_clustering(list(degree_centrality.values()), G, 'degree_centrality')
betweenness_kmeans, betweenness_optimal_k, betweenness_sse = perform_kmeans_clustering(list(betweenness_centrality.values()), G, 'betweenness_centrality')
closeness_kmeans, closeness_optimal_k, closeness_sse = perform_kmeans_clustering(list(closeness_centrality.values()), G, 'closeness_centrality')

centrality_df = pd.DataFrame({
    'degree_centrality': list(degree_centrality.values()),
    'betweenness_centrality': list(betweenness_centrality.values()),
    'closeness_centrality': list(closeness_centrality.values())
})

# Standardize the features
scaler = StandardScaler()
centrality_df_scaled = scaler.fit_transform(centrality_df)

# Elbow method 
sse_combined = []
for k in k_range:
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(centrality_df_scaled)
    sse_combined.append(kmeans.inertia_)

kneedle_combined = KneeLocator(k_range, sse_combined, curve='convex', direction='decreasing')
optimal_k_combined = kneedle_combined.elbow

# K-means clustering with the optimal number of clusters
kmeans_combined = KMeans(n_clusters=optimal_k_combined, random_state=42)
kmeans_combined.fit(centrality_df_scaled)

# Compute the mean centrality values for each cluster
centrality_means_combined = []
for cluster_label in range(optimal_k_combined):
    cluster_indices = np.where(kmeans_combined.labels_ == cluster_label)[0]
    mean_centrality = centrality_df.iloc[cluster_indices].mean().mean()
    centrality_means_combined.append((cluster_label, mean_centrality))

# Sort clusters by their mean centrality values
sorted_clusters_combined = sorted(centrality_means_combined, key=lambda x: x[1])

cluster_mapping_combined = {old_label: chr(65 + i) for i, (old_label, _) in enumerate(sorted_clusters_combined)}

# Add sorted cluster labels as node attributes
for i, node in enumerate(G.nodes):
    old_cluster_label = kmeans_combined.labels_[i]
    G.nodes[node]['combined_centrality_cluster'] = cluster_mapping_combined[old_cluster_label]

# Save the graph with clusters
nx.write_graphml(G, output_file)
print(f"Graph saved to {output_file}")

# Plot elbow graphs for each centrality measure
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(k_range, degree_sse, marker='o')
plt.vlines(degree_optimal_k, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
plt.xlabel('Number of clusters')
plt.ylabel('Sum of squared distances')
plt.title('Degree Centrality Elbow Curve')

plt.subplot(2, 2, 2)
plt.plot(k_range, betweenness_sse, marker='o')
plt.vlines(betweenness_optimal_k, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
plt.xlabel('Number of clusters')
plt.ylabel('Sum of squared distances')
plt.title('Betweenness Centrality Elbow Curve')

plt.subplot(2, 2, 3)
plt.plot(k_range, closeness_sse, marker='o')
plt.vlines(closeness_optimal_k, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
plt.xlabel('Number of clusters')
plt.ylabel('Sum of squared distances')
plt.title('Closeness Centrality Elbow Curve')

plt.subplot(2, 2, 4)
plt.plot(k_range, sse_combined, marker='o')
plt.vlines(optimal_k_combined, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
plt.xlabel('Number of clusters')
plt.ylabel('Sum of squared distances')
plt.title('Combined Centrality Elbow Curve')

plt.tight_layout()
plt.savefig(elbow_plot_file)
plt.close()

print(f"Elbow plots saved to {elbow_plot_file}")
