import pandas as pd

file = "data/raw/GSE184537_DiscoverySet_normalized_raw_counts.xlsx"

xls = pd.ExcelFile(file)

print(xls.sheet_names)