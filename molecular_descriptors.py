#Retrieve molecular structures and calculate descriptors

import requests
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
from rdkit.Chem import rdMolDescriptors

#Load cleaned potency data
df = pd.read_csv("egfr_cleaned_potency_data.csv")

compound_ids = df["molecule_chembl_id"].dropna().unique().tolist()

print("Number of compounds:", len(compound_ids))

#Retrieve molecular structures from ChEMBL
url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"

params = {
    "molecule_chembl_id__in": ",".join(compound_ids),
    "limit": 100
}

response = requests.get(url, params=params, timeout=60)

print("\nChEMBL RESPONSE CHECK")
print("Status code:", response.status_code)
print("Content type:", response.headers.get("Content-Type"))
print("Response preview:", response.text[:300])

response.raise_for_status()

molecules = response.json()["molecules"]

print("\nStructures retrieved:", len(molecules))

#Extract SMILES and calculate molecular descriptors
results = []

for molecule in molecules:

    chembl_id = molecule.get("molecule_chembl_id")

    structures = molecule.get("molecule_structures")

    if not structures:
        continue

    smiles = structures.get("canonical_smiles")

    if not smiles:
        continue

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        continue

    results.append({
        "molecule_chembl_id": chembl_id,
        "SMILES": smiles,
        "Molecular_Weight": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "Rotatable_Bonds": Lipinski.NumRotatableBonds(mol)
    })

#Create descriptor DataFrame
descriptor_df = pd.DataFrame(results)

#Merge descriptors with potency data
final_df = df.merge(
    descriptor_df,
    on="molecule_chembl_id",
    how="left"
)

#Save final dataset
final_df.to_csv("egfr_potency_with_descriptors.csv", index=False)

#Display results
print("\nMOLECULAR DESCRIPTOR ANALYSIS")

print("Compounds with structures:", len(descriptor_df))
print("Final dataset dimensions:", final_df.shape)

print("\nDESCRIPTOR SUMMARY")

print(
    descriptor_df[
        ["Molecular_Weight", "LogP", "HBD", "HBA",
         "TPSA", "Rotatable_Bonds"]
    ].describe().round(2)
)

print("\nTOP 5 COMPOUNDS BY pIC50")

print(
    final_df.sort_values("pIC50", ascending=False)[
        ["molecule_chembl_id", "pIC50", "Molecular_Weight",
         "LogP", "HBD", "HBA"]
    ].head(5).to_string(index=False)
)

print("\nSaved: egfr_potency_with_descriptors.csv")

