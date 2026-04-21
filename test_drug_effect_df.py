import json
import numpy as np
from neo4j import GraphDatabase
from tqdm import tqdm

def run_exhaustive_validation():
    with open('compound_rows.json', 'r') as f:
        row_labels = json.load(f)

    with open('protein_cols.json', 'r') as f:
        col_labels = json.load(f)

    matrix = np.load('dense_adjacency_matrix.npy')

    assert len(row_labels) == matrix.shape[0]
    assert len(col_labels) == matrix.shape[1]

    with open('compund_connectivity.json', 'r') as file:
        connectivity_data = json.load(file)

    rel_mapping = {item['connectionType']: item['matrixValue'] for item in connectivity_data}
    valid_rels = list(rel_mapping.keys())

    URI = "bolt://localhost:7683"
    AUTH = ("neo4j", "password")
    driver = GraphDatabase.driver(URI, auth=AUTH)

    def get_expected_db_value(tx, comp_id, prot_id):
        query = """
        MATCH (c:Compound)-[r]->(p:Protein)
        WHERE (coalesce(c.name, elementId(c)) = $cid) 
          AND (coalesce(p.name, elementId(p)) = $pid)
          AND type(r) IN $valid_rels
        RETURN type(r) AS rel_type
        """
        result = tx.run(query, cid=comp_id, pid=prot_id, valid_rels=valid_rels)
        records = result.data()
        
        if not records:
            return 0
        
        return max([rel_mapping[rec['rel_type']] for rec in records])

    non_zero_indices = np.argwhere(matrix != 0)
    
    with driver.session() as session:
        for r, c in tqdm(non_zero_indices, desc="Validating Non-Zero Edges"):
            comp_id = row_labels[r]
            prot_id = col_labels[c]
            matrix_val = matrix[r, c]
            db_val = session.execute_read(get_expected_db_value, comp_id, prot_id)
            assert matrix_val == db_val, f"Mismatch at ({r},{c}) for {comp_id}->{prot_id}. Matrix: {matrix_val}, DB: {db_val}"

    zero_indices = np.argwhere(matrix == 0)
    
    with driver.session() as session:
        for r, c in tqdm(zero_indices, desc="Validating Zero Edges"):
            comp_id = row_labels[r]
            prot_id = col_labels[c]
            matrix_val = matrix[r, c]
            db_val = session.execute_read(get_expected_db_value, comp_id, prot_id)
            assert matrix_val == db_val, f"Mismatch at ({r},{c}) for {comp_id}->{prot_id}. Matrix: {matrix_val}, DB: {db_val}"

    driver.close()
    print("Exhaustive matrix validation passed successfully.")

if __name__ == "__main__":
    run_exhaustive_validation()
