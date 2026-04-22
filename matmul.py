import numpy as np
import json
import os
import sys

def resort(list, idc):
    new_list = []
    for idx in idc:
        new_list.append(list[idx])
    return new_list


if len(sys.argv) != 6:
    print("Usage: python matmul.py <patient_protein_matrix.npy> <drug_weights.npy> <protein_drug_matrix.npy> <compound_rows.json> <output_drug_scores.json>")
    sys.exit(1)
os.makedirs(os.path.expanduser('./matrices/drug_scores'), exist_ok=True)
patient_protein_matrix_path = sys.argv[1]
drug_weights_path = sys.argv[2]
protein_drug_matrix_path = sys.argv[3]
compound_rows_path = sys.argv[4]
output_drug_scores_path = sys.argv[5]

patient_protein_matrix = np.load(os.path.expanduser(patient_protein_matrix_path))
drug_weights = np.load(os.path.expanduser(drug_weights_path))
protein_drug_matrix = np.load(os.path.expanduser(protein_drug_matrix_path)).transpose()

print("Patient-Features Matrix Shape:", patient_protein_matrix.shape)
print("Drug Weights Shape:", drug_weights.shape)
print("Features Drug Matrix Shape:", protein_drug_matrix.shape)

weighted_protein_drug_matrix = protein_drug_matrix #* drug_weights

patient_drug_matrix = np.matmul(patient_protein_matrix, weighted_protein_drug_matrix)
print("Patient-Drug Matrix Shape:", patient_drug_matrix.shape)

sorted_idx_drugs = np.argsort(patient_drug_matrix, axis=1)[0]



with open(os.path.expanduser(compound_rows_path), 'r') as f:
    compound_rows = json.load(f)
print(sorted_idx_drugs)
drug_scores = dict(zip(resort(compound_rows, sorted_idx_drugs), list(map(int, patient_drug_matrix[0, sorted_idx_drugs]))))
with open(os.path.expanduser(output_drug_scores_path), 'w') as f:
    json.dump(drug_scores, f, indent=4)
for i, idx in enumerate(sorted_idx_drugs):
    print("Sorted Drug Indices for a Patient:", idx, "Corresponding Drugs:", compound_rows[idx], "Patient-Drug Score:", patient_drug_matrix[0, idx])
    
    if i == 10:
        break