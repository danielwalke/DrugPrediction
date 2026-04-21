import numpy as np
import json

with open("matrices/gene_drug/compound_rows.json", 'r') as f:
    compound_rows = json.load(f)
with open("matrices/gene_drug/gene_cols.json", 'r') as f:
    gene_cols = json.load(f)

with open('matrices/drug_weights/drug_weights_genes.npy', 'rb') as f:
    drug_weights = np.load(f)

protein_drug_matrix = np.load('matrices/gene_drug/dense_adjacency_matrix.npy').transpose()*drug_weights
drug_idx = compound_rows.index('(3-ENDO)-8-METHYL-8-AZABICYCLO[3.2.1]OCT-3-YL 1H-PYRROLO[2,3-B]PYRIDINE-3-CARBOXYLATE')
regulated_protein_idc = np.where(protein_drug_matrix[:,drug_idx] != 0)
client_protein_matrix = np.load('matrices/patient_gene/patient_gene_matrix.npy')
print("Index of searched feature:", drug_idx)
print("Drug weight for searched feature:", drug_weights[drug_idx])
for idx in regulated_protein_idc[0]:
    print("Regulated gene:", gene_cols[idx])
print(client_protein_matrix[:,regulated_protein_idc])
print(protein_drug_matrix[regulated_protein_idc, drug_idx])