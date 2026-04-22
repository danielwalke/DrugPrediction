"""
Compute per-compound degree metrics used by drug_weights_calculation.py.

The drug prior w_d = 1 / sqrt(out_degree(d) + 1) is a *specificity*
prior: it penalizes compounds that interact with many entities across
the whole graph, since such promiscuity correlates with off-target
effects (other proteins, genes, transcripts, side effects, drug-drug
interactions, ...).

Edges to :Disease are somewhat ambiguous: a drug that treats many
diseases is not necessarily "unspecific" in the off-target sense; it
just has broad therapeutic use. Set the env var EXCLUDE_DISEASE=1 to
subtract :Disease out-edges from `out_degree` (and likewise in/complete
degree). Default: off, matching the legacy behaviour.
"""
from __future__ import annotations

import json
import os
import sys

from neo4j import GraphDatabase


EXCLUDE_DISEASE = os.environ.get("EXCLUDE_DISEASE", "0").lower() in {"1", "true", "yes"}


def process_compound_metrics(uri: str, user: str, password: str,
                             input_file: str, output_file: str) -> None:
    with open(input_file) as f:
        compound_names = json.load(f)

    if EXCLUDE_DISEASE:
        query = """
        UNWIND $names AS comp_name
        MATCH (c:Compound {name: comp_name})
        WITH c, comp_name,
             apoc.node.degree(c, '<') AS in_degree_raw,
             apoc.node.degree(c, '>') AS out_degree_raw,
             apoc.node.degree(c)      AS complete_degree_raw,
             COUNT { (c)-[]->(:Disease) } AS out_to_disease,
             COUNT { (:Disease)-[]->(c) } AS in_from_disease
        RETURN comp_name,
               in_degree_raw      - in_from_disease                     AS in_degree,
               out_degree_raw     - out_to_disease                      AS out_degree,
               complete_degree_raw - out_to_disease - in_from_disease   AS complete_degree,
               coalesce(c.betweennessCentrality, 0.0) AS betweenness
        """
    else:
        query = """
        UNWIND $names AS comp_name
        MATCH (c:Compound {name: comp_name})
        RETURN comp_name,
               apoc.node.degree(c, '<') AS in_degree,
               apoc.node.degree(c, '>') AS out_degree,
               apoc.node.degree(c)      AS complete_degree,
               coalesce(c.betweennessCentrality, 0.0) AS betweenness
        """

    driver = GraphDatabase.driver(uri, auth=(user, password))
    results_dict = {}

    with driver.session() as session:
        for record in session.run(query, names=compound_names):
            results_dict[record["comp_name"]] = {
                "in_degree": record["in_degree"],
                "out_degree": record["out_degree"],
                "complete_degree": record["complete_degree"],
                "betweenness_centrality": record["betweenness"],
            }

    driver.close()

    with open(output_file, "w") as f:
        json.dump(results_dict, f, indent=4)


if __name__ == "__main__":
    os.makedirs(os.path.expanduser("./matrices/drug_weights/"), exist_ok=True)
    if len(sys.argv) != 3:
        print("Usage: python degree-calculation.py <input_file_path> <output_file_path>")
        sys.exit(1)
    input_file_path, output_file_path = sys.argv[1], sys.argv[2]
    print(
        f"EXCLUDE_DISEASE={EXCLUDE_DISEASE}. "
        f"Processing compound metrics: {input_file_path} -> {output_file_path}"
    )
    process_compound_metrics(
        "bolt://localhost:7683",
        "neo4j",
        "password",
        os.path.expanduser(input_file_path),
        os.path.expanduser(output_file_path),
    )
