from neo4j import GraphDatabase
import json
import numpy as np
import os

os.makedirs(os.path.expanduser('./matrices/patient_gene/'), exist_ok=True)
driver = GraphDatabase.driver("bolt://localhost:7683", auth=("neo4j", "password"))
with open(os.path.expanduser('./matrices/gene_drug/gene_cols.json'), 'r') as f:
    gene_rows = json.load(f)

patient_gene_matrix = []
with driver.session() as session:
    result = session.run("""
        MATCH (p:Patient)-[r]->(q:Gene)
        RETURN p.name AS patient_name, q.name AS gene_name, r.regulation AS weight
        """).to_df().dropna()
result['weight'] = result['weight'].replace({"Down": -1, "Up": 1})

for patient_name, df_group in result.groupby('patient_name'):
    patient_gene_vector = []
    gene_name_to_weight_dict = dict(zip(df_group['gene_name'], df_group['weight']))
    for gene_name in gene_rows:
        patient_gene_vector.append(gene_name_to_weight_dict.get(gene_name, 0))
    patient_gene_matrix.append(patient_gene_vector)

np.save(os.path.expanduser('./matrices/patient_gene/patient_gene_matrix.npy'), np.array(patient_gene_matrix))


    