import pandas as pd
import os

# Read a file with obtained from cytoscape
df = pd.read_csv('../Files/STRING network default node.csv')

# Create a df with ENSG and ENSP columns only
df_pairs = df[["query term", "@id"]]
df_pairs = df_pairs.rename(columns={'query term': 'ENSG', '@id': 'ENSP'})

# Format ENSP column
for gene_pos in range(len(df_pairs)):
    df_pairs["ENSP"][gene_pos] = df_pairs["ENSP"][gene_pos][14:]

# Read a file with full ensembls of genes (with versions)
df_with_versions = pd.read_csv('../Files/ensembl_protein_coding.csv')

# Create a column without versions for merging
df_with_versions["ENSG"] = df_with_versions['gene_name'].apply(lambda x: x[:15])

# Merge the dataframes
df_pairs = pd.merge(df_pairs, df_with_versions, on="ENSG")
df_pairs.drop(columns=["ENSG"], inplace=True)
df_pairs = df_pairs.rename(columns={'gene_name': 'ENSG'})


folder_path = '../Files/clusters'

# Open each cluster
for filename in os.listdir(folder_path):

    with open(os.path.join(folder_path, filename), 'r') as cluster_file:

        ENSGs_to_write = []

        # Iterate through each gene of the cluster
        for ENSP in cluster_file:

            # Fix ENSP format
            ENSP = ENSP.strip()

            # Find a position of ENSP in df_pairs
            pos = df_pairs.index[df_pairs['ENSP'] == ENSP].tolist()

            # Get an ENSG from this position in df_pairs
            ENSG = df_pairs["ENSG"][pos[0]]

            # Add ENSG to the list
            ENSGs_to_write.append(ENSG)

        # Format the ENSGs before writing them to the output file
        ENSGs_to_write_str = '\n'.join(ENSGs_to_write)

        # Write the ENSGs to the file
        with open(f"../Files/clusters_ENSG/ENSG_{filename}", "w") as file_to_write:
            file_to_write.write(ENSGs_to_write_str)
