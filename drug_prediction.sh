python matrix_construction/protein_drug/create_drug_effect_df.py
python matrix_construction/drug_weights/degree-calculation.py "./matrices/protein_drug/compound_rows.json" "./matrices/drug_weights/compound_metrics.json"
python matrix_construction/drug_weights/drug_weights_calculation.py "./matrices/drug_weights/compound_metrics.json" "./matrices/protein_drug/compound_rows.json" "./matrices/drug_weights/drug_weights.npy"
python matrix_construction/patient_protein/construct_protein_matrix.py
python matmul.py "./matrices/patient_protein/patient_protein_matrix.npy" "./matrices/drug_weights/drug_weights.npy" "./matrices/protein_drug/dense_adjacency_matrix.npy" "./matrices/protein_drug/compound_rows.json" "./matrices/drug_scores/protein_drugs.json"
python matrix_construction/gene_drug/create_drug_effect_df.py
python matrix_construction/drug_weights/degree-calculation.py "./matrices/gene_drug/compound_rows.json" "./matrices/drug_weights/compound_metrics_genes.json"
python matrix_construction/drug_weights/drug_weights_calculation.py "./matrices/drug_weights/compound_metrics_genes.json" "./matrices/gene_drug/compound_rows.json" "./matrices/drug_weights/drug_weights_genes.npy"
python matrix_construction/gene_patient/construct_gene_matrix.py
python matmul.py "./matrices/patient_gene/patient_gene_matrix.npy" "./matrices/drug_weights/drug_weights_genes.npy" "./matrices/gene_drug/dense_adjacency_matrix.npy" "./matrices/gene_drug/compound_rows.json" "./matrices/drug_scores/gene_drugs.json"
python combined_score_analysis.py