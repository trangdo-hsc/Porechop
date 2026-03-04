import pandas as pd

barcodes_df = pd.read_csv("ONT_barcodes_SQK-NBD114.96.csv", sep=",", header=0)

barcodes_df.apply(lambda row: 
                  print(f"Adapter('Barcode {int(row.name) + 1} (forward)',\n\tstart_sequence=('{row['component']}', '{row['forward_sequence']}'),\n\tend_sequence=('{row['component']}_rev', '{row['reverse_sequence']}')),"), axis=1)