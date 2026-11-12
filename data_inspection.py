#Load chEMBL data
import pandas as pd
df = pd.read_csv("chembl_egfr_activity_data.csv")

#Disply dataset dimentions
print("\nDATASET SIZE")
print("Rows:", df.shape[0])
print("Columns:", df.shape[1])

#Display column names
print("\nCOLUMN NAMES")
for column in df.columns:
    print(column)

#Dsplay important experimental information
important_columns = ['molecule_chembl_id', 'assay_chembl_id', 'standard_type', 'standard_relation', 'standard_value', 'standard_units', 'pchembl_value']
print("\nIMPORTANT EXPERIMENTAL INFORMATION")
print(df[important_columns].head(10).to_string(index=False))

#Check for missing values in the dataset
print("\nMISSING VALUES")
print(df[important_columns].isnull().sum())
