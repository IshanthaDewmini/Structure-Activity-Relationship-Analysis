#EGFR inhibitor SAR Analysis
import pandas as pd
#Load potency data with molecular descriptors
df = pd.read_csv("egfr_potency_with_descriptors.csv")
#sort compounds by potency
top_compounds = df.sort_values(by="pIC50", ascending=False).head(5)
#Display important information for the top 5 compounds
print("\nTOP 5 COMPOUNDS BY pIC50")
print(
    top_compounds[
        ["molecule_chembl_id", "pIC50", "Molecular_Weight",
         "LogP", "HBD", "HBA", "TPSA", "SMILES"]
    ].to_string(index=False)
)
print("\nNumber of compounds selected:", len(top_compounds))


#Visualise the five most potent compounds
from rdkit import Chem
from rdkit.Chem import Draw
molecules = []
legends = []
for _, row in top_compounds.iterrows():

    mol = Chem.MolFromSmiles(row["SMILES"])

    if mol is not None:

        molecules.append(mol)

        legends.append(
            f'{row["molecule_chembl_id"]}\n'
            f'pIC50 = {row["pIC50"]:.2f}'
        )
#Generate molecular structure image
image = Draw.MolsToGridImage(
    molecules,
    molsPerRow=2,
    subImgSize=(450, 350),
    legends=legends
)
#Save image
image.save("top5_egfr_structures.png")

print("\nMolecular structures saved:")
print("top5_egfr_structures.png")
