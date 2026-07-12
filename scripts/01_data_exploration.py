import pandas as pd

# Load Excel workbook
file = "data/raw/GSE184537_DiscoverySet_normalized_raw_counts.xlsx"

xls = pd.ExcelFile(file)

print("\nAvailable sheets:")
print(xls.sheet_names)

# Read expression matrix
df = pd.read_excel(file, sheet_name=xls.sheet_names[0])

print("\nDataset shape:")
print(df.shape)

print("\nFirst five rows:")
print(df.head())

print("\nColumn names:")
print(df.columns.tolist())

print("\nMissing values:")
print(df.isna().sum().sum())
