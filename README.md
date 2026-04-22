# DrugPrediction

## Setup
1. Download latest pharmebidump here:
https://pharmebi.net/#/download
2. Create database dump:
sh pharmebinet/create_db.sh
3. Start database
sh pharmebinet/start_db.sh

python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

## Protein based prediction
1. Construct protein drug matrix (PxD) matrix 
python matrix_construction/protein_drug/create_drug_effect_df.py

2. Create Drug degree json
python matrix_construction/drug_weights/degree-calculation.py  "./matrices/protein_drug/compound_rows.json" "./matrices/drug_weights/compound_metrics.json"

3. Based on these degree we calculate the drug weight
python matrix_construction/drug_weights/drug_weights_calculation.py  './matrices/drug_weights/compound_metrics.json' './matrices/protein_drug/compound_rows.json' './matrices/drug_weights/drug_weights.npy'

4. Construct the protein to patient matrix 
python matrix_construction/patient_protein/construct_protein_matrix.py 

5. Use matrix multiplication get the best drugs (mathemagically)
python matmul.py './matrices/patient_protein/patient_protein_matrix.npy' './matrices/drug_weights/drug_weights.npy' './matrices/protein_drug/dense_adjacency_matrix.npy' './matrices/protein_drug/compound_rows.json' './matrices/drug_scores/protein_drugs.json'


## Gene based predcition
1. Construct gene drug matrix (PxD) matrix
python matrix_construction/gene_drug/create_drug_effect_df.py

2. Create Drug degree matrix
python matrix_construction/drug_weights/degree-calculation.py  "./matrices/gene_drug/compound_rows.json" "./matrices/drug_weights/compound_metrics_genes.json"

3. Based on these degree we calculate the drug weight
python matrix_construction/drug_weights/drug_weights_calculation.py  './matrices/drug_weights/compound_metrics_genes.json' './matrices/gene_drug/compound_rows.json' './matrices/drug_weights/drug_weights_genes.npy'

4. Construct the gene to patient matrix 
python matrix_construction/gene_patient/construct_gene_matrix.py 

5. Use matrix multiplication get the best drugs (mathemagically)
python matmul.py './matrices/patient_gene/patient_gene_matrix.npy' './matrices/drug_weights/drug_weights_genes.npy' './matrices/gene_drug/dense_adjacency_matrix.npy' './matrices/gene_drug/compound_rows.json' './matrices/drug_scores/gene_drugs.json'

## Combine analyiss (only available after running the previous commands)
python combined_score_analysis.py 

Outputs a json in 
"./matrics/drug_scores/combined_drugs.json"