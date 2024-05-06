import os

folder_path = '../Files/clusters'
file_gene_counts = {}

# Open each cluster
for filename in os.listdir(folder_path):
    with open(os.path.join(folder_path, filename), 'r') as file:
        line_count = sum(1 for line in file)
        # Store the file name and its line count in the dictionary
        file_gene_counts[filename] = line_count

# Sort the dictionary by gene count in descending order
sorted_files = sorted(file_gene_counts.items(), key=lambda x: x[1], reverse=True)

# Get the top 10 files with the most genes
top_10_files = sorted_files[:10]

print(top_10_files)