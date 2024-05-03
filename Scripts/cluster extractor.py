clusters = []
count = 0

# Open Clusters.txt for reading
with open('../Files/Clusters.txt', 'r') as file:
    # Iterate through each line to find the first line that contains information about cluster
    for line in file:
        if "Cluster	Score" in line:
            break

    unique_ids = []
    # Iterate through each line that contains information about cluster
    for line in file:
        count +=1   # For file names

        # Split the line by ', ' into parts
        parts = line.split(', ')
        ids = []

        # Iterate through each part in the line to get IDs
        for part in parts:
            # Add all genes to the list that will be saved to the file of cluster
            if '9606' in part:
                gene_to_add = part.split('.')[1]
                ids.append(gene_to_add)

                # Add all unique genes to the list that will be saved to the file with all genes
                if gene_to_add not in unique_ids:
                    unique_ids.append(gene_to_add)

        # Format the IDs before writing them to the output file
        ids_str = '\n'.join(ids)

        # Write the IDs to the file
        with open(f"../Files/clusters/cluster{count}.txt", "w") as file_to_write:
            file_to_write.write(ids_str)


# Remove '\n' characters after gene name at the end of the cluster
unique_ids_without_newline = [id.rstrip('\n') for id in unique_ids]

# Format the IDs before writing them to the output file
unique_ids_without_newline_str = '\n'.join(unique_ids_without_newline)

# Write all the IDs to one file
with open(f"../Files/DEGs_from_clusters.txt", "w") as file_to_write:
    file_to_write.write(unique_ids_without_newline_str)