#Identification of largest assay
import pandas as pd
# Load the ChEMBL dataset
df = pd.read_csv("chembl_egfr_activity_data.csv")
# Select the largest assay
assay_id = "CHEMBL674643"

assay_data = df[df["assay_chembl_id"] == assay_id].copy()

print("\nASSAY INFORMATION")
print("Assay ID:", assay_id)
print("Number of records:", len(assay_data))

print("\nPUBLICATION INFORMATION")

publication_columns = [
    "document_chembl_id",
    "document_journal",
    "document_year",
    "assay_description"
]
print(
    assay_data[publication_columns]
    .drop_duplicates()
    .to_string(index=False)
)
print("\nCOMPOUNDS TESTED")

compound_columns = [
    "molecule_chembl_id",
    "standard_value",
    "standard_units",
    "pchembl_value"
]
print(
    assay_data[compound_columns]
    .to_string(index=False)
)

#Find the manuscript details for the assay
import requests
document_id = "chembl1124610"
url = f"https://www.ebi.ac.uk/chembl/api/data/document/{document_id}.json"
response = requests.get(url, timeout=30)
response.raise_for_status()
paper = response.json()
print("\nORIGINAL RESEARCH PAPER DETAILS")
print("Title:", paper.get("title"))
print("Authors:", paper.get("authors"))
print("Journal:", paper.get("journal"))
print("Year:", paper.get("year"))
print("DOI:", paper.get("doi"))
print("PubMed ID:", paper.get("pubmed_id"))
