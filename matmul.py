import numpy as np
import json
import os
import sys

if len(sys.argv) != 5:
    print("Usage: python matmul.py <patient_protein_matrix.npy> <drug_weights.npy> <protein_drug_matrix.npy> <compound_rows.json>")
    sys.exit(1)

patient_protein_matrix_path = sys.argv[1]
drug_weights_path = sys.argv[2]
protein_drug_matrix_path = sys.argv[3]
compound_rows_path = sys.argv[4]

patient_protein_matrix = np.load(os.path.expanduser(patient_protein_matrix_path))
drug_weights = np.load(os.path.expanduser(drug_weights_path))
protein_drug_matrix = np.load(os.path.expanduser(protein_drug_matrix_path)).transpose()

print("Patient-Features Matrix Shape:", patient_protein_matrix.shape)
print("Drug Weights Shape:", drug_weights.shape)
print("Features Drug Matrix Shape:", protein_drug_matrix.shape)



weighted_protein_drug_matrix = protein_drug_matrix * drug_weights

patient_drug_matrix = np.matmul(patient_protein_matrix, weighted_protein_drug_matrix)
print("Patient-Drug Matrix Shape:", patient_drug_matrix.shape)

sorted_idx_drugs = np.argsort(patient_drug_matrix, axis=1)[0]
print(patient_drug_matrix)
print(sorted_idx_drugs)
print(patient_drug_matrix[:, sorted_idx_drugs])



with open(os.path.expanduser(compound_rows_path), 'r') as f:
    compound_rows = json.load(f)
print(sorted_idx_drugs)
for i, idx in enumerate(sorted_idx_drugs):
    print("Sorted Drug Indices for a Patient:", idx, "Corresponding Drugs:", compound_rows[idx], "Patient-Drug Score:", patient_drug_matrix[0, idx])
    if i == 10:
        break