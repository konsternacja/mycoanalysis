from biom import load_table
from sys import argv

# Load the first table
merged_table = load_table(argv[1])

# Merge with the remaining tables
for table_path in argv[2:-1]:
    table = load_table(table_path)
    merged_table = merged_table.merge(table)

# Write the merged table to a file
with open(argv[-1], 'w') as f:
    f.write(merged_table.to_json("Merged table"))