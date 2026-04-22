"""
Compute per-compound degree metrics used by drug_weights_calculation.py.

The drug prior w_d = 1 / sqrt(out_degree(d) + 1) is a *specificity*
prior: compounds that interact with many other biological entities are
down-weighted because such promiscuity correlates with off-target
effects (drug-drug interactions, non-target protein/gene binding, ...).

Which labels should count towards that degree is a modelling choice,
controlled by the DEGREE_SCOPE environment variable:

    DEGREE_SCOPE=targets  (default)
        Count only edges to :Compound, :Protein and :Gene nodes. This
        is the right choice if you want the specificity prior to
        reflect *interaction* promiscuity: drug-drug interactions,
        protein binding, gene regulation. Edges to :Disease,
        :SideEffect, :Pathway etc. are excluded because e.g. "treats
        many diseases" does not imply "hits many targets".

    DEGREE_SCOPE=all
        Count edges to every connected node, matching the legacy
        apoc.node.degree behaviour.
"""
from __future__ import annotations

import json
import os
import sys

from neo4j import GraphDatabase


DEGREE_SCOPE = os.environ.get("DEGREE_SCOPE", "targets").lower()
if DEGREE_SCOPE not in {"targets", "all"}:
    raise ValueError(
        f"DEGREE_SCOPE must be one of targets|all, got {DEGREE_SCOPE!r}"
    )


QUERY_TARGETS = """
UNWIND $names AS comp_name
MATCH (c:Compound {name: comp_name})
WITH c, comp_name,
     COUNT { (c)-[]->(:Compound) }
   + COUNT { (c)-[]->(:Protein)  }
   + COUNT { (c)-[]->(:Gene)     }                  AS out_degree,
     COUNT { (:Compound)-[]->(c) }
   + COUNT { (:Protein)-[]->(c)  }
   + COUNT { (:Gene)-[]->(c)     }                  AS in_degree
RETURN comp_name,
       in_degree,
       out_degree,
       in_degree + out_degree                        AS complete_degree,
       coalesce(c.betweennessCentrality, 0.0)        AS betweenness
"""

QUERY_ALL = """
UNWIND $names AS comp_name
MATCH (c:Compound {name: comp_name})
RETURN comp_name,
       apoc.node.degree(c, '<') AS in_degree,
       apoc.node.degree(c, '>') AS out_degree,
       apoc.node.degree(c)      AS complete_degree,
       coalesce(c.betweennessCentrality, 0.0) AS betweenness
"""


def process_compound_metrics(uri: str, user: str, password: str,
                             input_file: str, output_file: str) -> None:
    with open(input_file) as f:
        compound_names = json.load(f)

    query = QUERY_TARGETS if DEGREE_SCOPE == "targets" else QUERY_ALL

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
        f"DEGREE_SCOPE={DEGREE_SCOPE}. "
        f"Processing compound metrics: {input_file_path} -> {output_file_path}"
    )
    process_compound_metrics(
        "bolt://localhost:7683",
        "neo4j",
        "password",
        os.path.expanduser(input_file_path),
        os.path.expanduser(output_file_path),
    )
