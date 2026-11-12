print("ChEMBL Dose-Response Analysis")
import requests
import pandas as pd
url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
params = {
    "limit": 100,
    "target_chembl_id" : "CHEMBL203",
    "standard_type": "IC50",
}
#Retrieve Data from chEMBL
response = requests.get(url, params=params, timeout=30)
response.raise_for_status()
data = response.json()
#Convert activity data to a pandas DataFrame
activities = data['activities']
df = pd.DataFrame(activities)
#Display information about the DataFrame
print("EGFR activity data retrieved from ChEMBL")
print("Number of records:", len(df))
print("\nFirst five redords:")
print(df.head())
#Save the original data to a CSV file
df.to_csv("chembl_egfr_activity_data.csv", index=False)
print("\nData saved to 'chembl_egfr_activity_data.csv'")