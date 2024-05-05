clusters = []
count = 0

# Open Clusters.txt for reading
with open('../Files/Clusters.txt', 'r') as file:
    # Iterate through each line to find the first line that contains information about cluster
    for line in file:
        if "Cluster	Score" in line:
            break

    # Iterate through each line that contains information about cluster
    for line in file:
        count +=1   # For file names

        # Split the line by ', ' into parts
        parts = line.split(', ')
        ids = []

        # Iterate through each part in the line to get IDs
        for part in parts:
            if '9606' in part:
                ids.append(part.split('.')[1])

        # Format the IDs before writing them to the output file
        ids_str = '\n'.join(ids)

        # Write the IDs to the file
        with open(f"../Files/clusters/cluster{count}.txt", "w") as file_to_write:
            file_to_write.write(ids_str)

