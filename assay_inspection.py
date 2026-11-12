#Identification of unique assays

import pandas as pd
#Load chEMBL data
df = pd.read_csv("chembl_egfr_activity_data.csv")
#Display unique assay identifiers
print("\nUNIQUE ASSAY IDENTIFIERS")
print("Number of unique assay identifiers:", df['assay_chembl_id'].nunique())

#Summarise Each assay by counting the number of records associated with it
assay_summary = df.groupby("assay_chembl_id").agg(
    number_of_records=("activity_id", "count"),
    assay_type=("assay_type", "first"),
    assay_description=("assay_description", "first"),
    target_name=("target_pref_name", "first"),
    organism=("target_organism", "first"),
    document_year=("document_year", "first")
).reset_index()
print("\nASSAY SUMMARY")
print(assay_summary.to_string(index=False))

#Save summary
assay_summary.to_csv(
    "assay_summary.csv",
    index=False
)
print("\nAssay summary saved successfully!")