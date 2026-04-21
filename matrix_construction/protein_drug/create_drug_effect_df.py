import json
import pandas as pd
import numpy as np
from neo4j import GraphDatabase
from scipy.sparse import csr_matrix, save_npz
import os

os.makedirs(os.path.expanduser('./matrices/protein_drug/'), exist_ok=True)
with open(os.path.expanduser('./config/compound_connectivity.json'), 'r') as file:
    connectivity_data = json.load(file)

rel_mapping = {item['connectionType']: item['matrixValue'] for item in connectivity_data}
rel_types = list(rel_mapping.keys())

URI = "bolt://localhost:7683"
AUTH = ("neo4j", "password")

driver = GraphDatabase.driver(URI, auth=AUTH)

query = """
MATCH (c:Compound)-[r]->(p:Protein)
WHERE type(r) IN $rel_types
RETURN coalesce(c.name, elementId(c)) AS compound_id, type(r) AS connection_type, coalesce(p.name, elementId(p)) AS protein_id
ORDER BY elementId(r)
SKIP $skip LIMIT $limit
"""

batch_size = 10000
skip = 0
all_records = []

with driver.session() as session:
    while True:
        result = session.run(query, rel_types=rel_types, skip=skip, limit=batch_size)
        records = result.data()
        
        if not records:
            break
            
        all_records.extend(records)
        skip += batch_size

driver.close()

df = pd.DataFrame(all_records)

if not df.empty:
    df['matrix_value'] = df['connection_type'].map(rel_mapping)

    matrix_df = df.pivot_table(
        index='compound_id', 
        columns='protein_id', 
        values='matrix_value', 
        fill_value=0,
        aggfunc='max'
    )

    matrix_df['row_sum'] = matrix_df.sum(axis=1)
    matrix_df.loc['col_sum'] = matrix_df.sum(axis=0)

    pure_matrix = matrix_df.drop(index='col_sum', columns='row_sum')

    row_labels = pure_matrix.index.tolist()
    col_labels = pure_matrix.columns.tolist()

    with open(os.path.expanduser('./matrices/protein_drug/compound_rows.json'), 'w') as f:
        json.dump(row_labels, f)

    with open(os.path.expanduser('./matrices/protein_drug/protein_cols.json'), 'w') as f:
        json.dump(col_labels, f)

    numpy_matrix = pure_matrix.to_numpy(dtype=np.int8)
    np.save(os.path.expanduser('./matrices/protein_drug/dense_adjacency_matrix.npy'), numpy_matrix)

    sparse_matrix = csr_matrix(numpy_matrix)
    save_npz(os.path.expanduser('./matrices/protein_drug/sparse_adjacency_matrix.npz'), sparse_matrix)

    cols_to_print = list(matrix_df.columns[:10]) + ['row_sum']
    rows_to_print = list(matrix_df.index[:6]) + ['col_sum']
    
    print("--- Matrix Snapshot (6x10 with Sums) ---")
    print(matrix_df.loc[rows_to_print, cols_to_print])
    print("\n")

    print("--- First 10 Examples of Upregulation / Agonism (+1) ---")
    print(df[df['matrix_value'] == 1][['compound_id', 'protein_id', 'connection_type']].head(10))
    print("\n")

    print("--- First 10 Examples of Downregulation / Inhibition (-1) ---")
    print(df[df['matrix_value'] == -1][['compound_id', 'protein_id', 'connection_type']].head(10))

else:
    print("Empty DataFrame. No relationships found.")
