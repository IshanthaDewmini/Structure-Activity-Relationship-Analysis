#Clean the potency data and calculate pIC50
import pandas as pd
import numpy as np
df = pd.read_csv("chembl_egfr_activity_data.csv")
print("Original dataset:", df.shape)
#Select one assay and IC50 data
assay_df = df[
    (df["assay_chembl_id"] == "CHEMBL674643") &
    (df["standard_type"] == "IC50")
].copy()
print("Records from selected assay:", len(assay_df))

#Convert activity values to numeric
assay_df["standard_value"] = pd.to_numeric(assay_df["standard_value"], errors="coerce")
#Keep valid measurements
clean_df = assay_df[
    (assay_df["standard_units"] == "nM") &
    (assay_df["standard_relation"] == "=") &
    (assay_df["standard_value"].notna()) &
    (assay_df["standard_value"] > 0)
].copy()

print("Valid IC50 measurements:", len(clean_df))

#Calculate pIC50 values
clean_df["pIC50"] = -np.log10(
    clean_df["standard_value"] * 1e-9
)

#Check duplicate compound measurements
duplicates = clean_df.duplicated(subset=["molecule_chembl_id", "assay_chembl_id"]).sum()
print("Duplicate compound measurements:", duplicates)

#Display most potent compounds
clean_df = clean_df.sort_values("pIC50", ascending=False)
print("\nMOST POTENT COMPOUNDS")
print(clean_df[["molecule_chembl_id", "standard_value", "pIC50", "standard_units"]].head(10).to_string(index=False))

#Summary statistics of pIC50 values
print("\nPOTENCY SUMMARY")
print(clean_df["pIC50"].describe())

#Save cleaned data
clean_df.to_csv("egfr_cleaned_potency_data.csv", index=False)
print("\nCleaned dataset saved successfully.")
                