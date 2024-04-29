
with open("../Files/ensembl_protein_coding.csv") as file:
    ensembles = []

    for line in file:
        ensembles.append(line[1:16])

ensembles.remove('gene_name"\n')


# Format the IDs before writing them to the output file
ids_str = '\n'.join(ensembles)

# Write the IDs to the file
with open("../Files/ensembl_without_versions.csv", "w") as file_to_write:
    file_to_write.write(ids_str)


# # THERE IS NO DUPLICATES
# print(len(ensembles))
# ensembles_unique = []
# for i in ensembles:
#     if i not in ensembles_unique:
#         ensembles_unique.append(i)
# print(len(ensembles_unique))