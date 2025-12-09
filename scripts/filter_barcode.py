import pandas as pd
import sys

# Define DataFrame to filter from argument path, read csv
df_path = sys.argv[1]
cell_barcodes = pd.read_csv(df_path)

# Define filtering parameters
cell_type = sys.argv[2]
sample = sys.argv[3]
print(cell_barcodes[
    (cell_barcodes['cell_type'] == cell_type) & 
    (cell_barcodes['sample'] == sample)].shape)
# Where to save the output
output_path = sys.argv[4]

barcodes = cell_barcodes[
    (cell_barcodes['cell_type'] == cell_type) & 
    (cell_barcodes['sample'] == sample)]['barcode'].to_list()
print(barcodes)

file = open(output_path, 'w')
for x in ['CB:Z:' + x + '\n' for x in barcodes]:
    file.write(x)
file.close()