from neo4j import GraphDatabase
import json
import numpy as np
import os

os.makedirs(os.path.expanduser('./matrices/patient_protein/'), exist_ok=True)
driver = GraphDatabase.driver("bolt://localhost:7683", auth=("neo4j", "password"))
with open(os.path.expanduser('./matrices/protein_drug/protein_cols.json'), 'r') as f:
    protein_rows = json.load(f)

patient_protein_matrix = []
with driver.session() as session:
    result = session.run("""
        MATCH (p:Patient)-[r]->(q:Protein)
        RETURN p.name AS patient_name, q.name AS protein_name, r.regulation AS weight
        """).to_df().dropna()
result['weight'] = result['weight'].replace({"Down": -1, "Up": 1})

for patient_name, df_group in result.groupby('patient_name'):
    patient_protein_vector = []
    protein_name_to_weight_dict = dict(zip(df_group['protein_name'], df_group['weight']))
    for protein_name in protein_rows:
        patient_protein_vector.append(protein_name_to_weight_dict.get(protein_name, 0))
    patient_protein_matrix.append(patient_protein_vector)

np.save(os.path.expanduser('./matrices/patient_protein/patient_protein_matrix.npy'), np.array(patient_protein_matrix))



    