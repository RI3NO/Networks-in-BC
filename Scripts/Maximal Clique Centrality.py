import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def read_sif_file(file_path):
    """Read the .sif file and create the graph"""
    G = nx.Graph()
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                G.add_edge(parts[0], parts[2])
            elif len(parts) == 1:
                G.add_node(parts[0])
    return G

def calculate_maximal_clique_centrality(graph):
    """Calculate Maximal Clique Centrality for each node of the graph"""
    cliques = list(nx.find_cliques(graph))
    clique_count = {node: 0 for node in graph.nodes()}

    for clique in cliques:
        for node in clique:
            clique_count[node] += 1

    return clique_count

def normalize_centrality(clique_count):
    """Normalize Maximal Clique Centrality"""
    max_centrality = max(clique_count.values()) if clique_count else 0
    normalized_centrality = {node: count / max_centrality if max_centrality > 0 else 0 for node, count in clique_count.items()}
    return normalized_centrality


# Create the graph from the .sif file
G = read_sif_file("../Files/SHIT network.sif")
# G = read_sif_file("../Files/STRING results.sif")

# Finding all maximal cliques
cliques = list(nx.find_cliques(G))
largest_clique = max(cliques, key=len)

print("All cliques:", cliques)
print("Largest clique:", largest_clique)


# Calculate Maximal Clique Centrality
clique_count = calculate_maximal_clique_centrality(G)
normalized_centrality = normalize_centrality(clique_count)

# Transform to Pandas Series
clique_centrality_series = pd.Series(normalized_centrality, name='MCC')
clique_centrality_series_normalized = pd.Series(clique_count, name='MCC')

print(clique_centrality_series.sort_values(ascending=False))
print(clique_centrality_series_normalized.sort_values(ascending=False))


# Merge the two Series into a DataFrame
df = pd.DataFrame({
    'MCC': clique_centrality_series_normalized,
    'Normalized MCC': clique_centrality_series
})

# Sort the DataFrame by 'MCC' in descending order
df_sorted = df.sort_values(by='MCC', ascending=False)

# Save the sorted DataFrame to a CSV file
df_sorted.to_csv('../Files/sorted_clique_centrality.csv')


# # Graph drawing
# pos = nx.spring_layout(G)  # positions for all nodes
#
# # Draw all nodes and edges
# nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=500)
# nx.draw_networkx_edges(G, pos, alpha=0.5)
#
# # Highlight the nodes in the largest clique
# nx.draw_networkx_nodes(G, pos, nodelist=largest_clique, node_color='orange', node_size=500)
#
# # Labels for the nodes
# # nx.draw_networkx_labels(G, pos, font_size=16, font_family='sans-serif')
#
# plt.title('Graph with Largest Clique Highlighted')
# plt.axis('off')  # Turn off the axis
# plt.show()








# # TEST TEST TEST
# # Create a graph
# G = nx.Graph()
#
# # Add nodes
# G.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G'])
#
# # Add edges
# edges = [
#     ('A', 'B'), ('A', 'C'), ('A', 'D'), ('A', 'F'), ('A', 'G'),
#     ('B', 'C'), ('B', 'E'), ('B', 'D'),
#     ('C', 'D'), ('C', 'F'), ('C', 'G'), ('C', 'E'),
#     ('D', 'E'), ('D', 'F'), ('D', 'G'),
#     ('E', 'F'), ('E', 'G'),
#     ('F', 'G'),
# ]
# G.add_edges_from(edges)
#
# TESTTESTTEST = pd.read_csv('../Files/some random shit.csv')
# print(TESTTESTTEST.sort_index())
#
#
# # Finding all maximal cliques
# cliques = list(nx.find_cliques(G))
# largest_clique = max(cliques, key=len)
#
# print("All cliques:", cliques)
# print("Largest clique:", largest_clique)
#
#
# # Calculate Maximal Clique Centrality
# clique_count = calculate_maximal_clique_centrality(G)
# normalized_centrality = normalize_centrality(clique_count)
#
# # Transform to Pandas Series
# clique_centrality_series = pd.Series(normalized_centrality, name='MCC')
# clique_centrality_series_normalized = pd.Series(clique_count, name='MCC')
#
# print(clique_centrality_series.sort_index())
# print(clique_centrality_series_normalized.sort_index())