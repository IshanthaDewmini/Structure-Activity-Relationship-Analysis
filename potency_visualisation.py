#EGFR Inhibitor Potency and Molecular Property Analysis
import pandas as pd
import matplotlib.pyplot as plt

#Load potency data with molecular descriptors
df = pd.read_csv("egfr_potency_with_descriptors.csv")
print("Dataset dimensions:", df.shape)
print("\nAvailable columns:")
print(df.columns.tolist())

#Define Molecular descriptors to visualize
descriptors = [
    "Molecular_Weight",
    "LogP",
    "HBD",
    "HBA",
    "TPSA",
    "Rotatable_Bonds"
]

#Check for missing values
print("\nMissing values:")
print(df[["pIC50"] + descriptors].isnull().sum())

#Plot pIC50 distribution
plt.figure(figsize=(8, 5))
plt.hist(df["pIC50"].dropna(), bins=30, color="skyblue", edgecolor="black")
plt.title("Distribution of EGFR Inhibitor Potency (pIC50)")
plt.xlabel("pIC50")
plt.ylabel("Number of Compounds-Frequency")
plt.tight_layout()
plt.savefig("EGFR_Inhibitor_pIC50_Distribution.png", dpi=300)
plt.show()

#Plot molecular descriptor distributions
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
for ax, descriptor in zip(axes.flatten(), descriptors):
    ax.scatter(df[descriptor], df["pIC50"])
    ax.set_xlabel(descriptor)
    ax.set_ylabel("pIC50")
plt.tight_layout()
plt.savefig("EGFR_Inhibitor_Descriptor_Plots.png", dpi=300)
plt.show()

#Calculate pearson correlation coefficients between pIC50 and molecular descriptors
correlations = df[["pIC50"] + descriptors].corr(
    method="pearson"
)
print("\nPEARSON CORRELATION WITH pIC50")
print(
    correlations["pIC50"]
    .drop("pIC50")
    .sort_values(key=abs, ascending=False)
    .round(3)
)
print("\nAnalysis completed successfully!")