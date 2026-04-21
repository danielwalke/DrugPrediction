import numpy as np
import json
import os

patient_protein_matrix = np.load(os.path.expanduser('./matrices/patient_protein/patient_protein_matrix.npy'))
drug_weights = np.load(os.path.expanduser('./matrices/drug_weights/drug_weights.npy'))
protein_drug_matrix = np.load(os.path.expanduser('./matrices/protein_drug/dense_adjacency_matrix.npy')).transpose()

weighted_protein_drug_matrix = protein_drug_matrix * drug_weights

patient_drug_matrix = np.matmul(patient_protein_matrix, weighted_protein_drug_matrix)
print("Patient-Protein Matrix Shape:", patient_protein_matrix.shape)
print("Drug Weights Shape:", drug_weights.shape)
print("Protein-Drug Matrix Shape:", protein_drug_matrix.shape)
print("Patient-Drug Matrix Shape:", patient_drug_matrix.shape)

sorted_idx_drugs = np.argsort(patient_drug_matrix, axis=1)[0]
print(patient_drug_matrix)
print(sorted_idx_drugs)
print(patient_drug_matrix[:, sorted_idx_drugs])



with open(os.path.expanduser('./matrices/protein_drug/compound_rows.json'), 'r') as f:
    compound_rows = json.load(f)
print(sorted_idx_drugs)
for i, idx in enumerate(sorted_idx_drugs):
    print("Sorted Drug Indices for a Patient:", idx, "Corresponding Drugs:", compound_rows[idx], "Patient-Drug Score:", patient_drug_matrix[0, idx])
    if i == 10:
        break