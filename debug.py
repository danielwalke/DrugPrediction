import numpy as np
import json

with open("compound_rows.json", 'r') as f:
    compound_rows = json.load(f)
with open("protein_cols.json", 'r') as f:
    protein_cols = json.load(f)

with open('drug_weights.npy', 'rb') as f:
    drug_weights = np.load(f)

protein_drug_matrix = np.load('dense_adjacency_matrix.npy').transpose()*drug_weights
drug_idx = compound_rows.index('CRA-028129')
regulated_protein_idc = np.where(protein_drug_matrix[:,drug_idx] != 0)
client_protein_matrix = np.load('patient_protein_matrix.npy')
print("Index of CRA-028129:", drug_idx)
print("Drug weight for CRA-028129:", drug_weights[drug_idx])
for idx in regulated_protein_idc[0]:
    print("Regulated protein:", protein_cols[idx])
print(client_protein_matrix[:,regulated_protein_idc])
print(protein_drug_matrix[regulated_protein_idc, drug_idx])