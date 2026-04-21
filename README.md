1. Construct protein drug matrix (PxD) matrix 
python matrix_construction/protein_drug/create_drug_effect_df.py

2. Create Drug degree json
python matrix_construction/drug_weights/degree-calculation.py 

3. Based on these degree we calculate the drug weight
python matrix_construction/drug_weights/drug_weights_calculation.py 

4. Construct the protein to patient matrix 
python matrix_construction/patient_protein/construct_protein_matrix.py 

5. Use matrix multiplication get the best drugs (mathemagically)
python matmul.py# DrugPrediction
