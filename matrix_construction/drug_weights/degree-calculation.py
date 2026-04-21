import json
from neo4j import GraphDatabase
import os

def process_compound_metrics(uri, user, password, input_file, output_file):
    with open(input_file, 'r') as f:
        compound_names = json.load(f)

    driver = GraphDatabase.driver(uri, auth=(user, password))
    
    query = """
    UNWIND $names AS comp_name
    MATCH (c:Compound {name: comp_name})
    RETURN comp_name,
           apoc.node.degree(c, '<') AS in_degree,
           apoc.node.degree(c, '>') AS out_degree,
           apoc.node.degree(c) AS complete_degree,
           coalesce(c.betweennessCentrality, 0.0) AS betweenness
    """

    results_dict = {}

    with driver.session() as session:
        result = session.run(query, names=compound_names)
        for record in result:
            results_dict[record["comp_name"]] = {
                "in_degree": record["in_degree"],
                "out_degree": record["out_degree"],
                "complete_degree": record["complete_degree"],
                "betweenness_centrality": record["betweenness"]
            }

    driver.close()

    with open(output_file, 'w') as f:
        json.dump(results_dict, f, indent=4)

if __name__ == "__main__":
    os.makedirs(os.path.expanduser('./matrices/drug_weights/'), exist_ok=True)
    process_compound_metrics(
        "bolt://localhost:7683",
        "neo4j",
        "password",
        os.path.expanduser("./matrices/protein_drug/compound_rows.json"),
        os.path.expanduser("./matrices/drug_weights/compound_metrics.json")
    )
