"""
Build the patient x gene regulation matrix.

See matrix_construction/patient_protein/construct_protein_matrix.py for
the weighting rationale; this file is the gene analogue.
"""
from __future__ import annotations

import json
import os

import numpy as np
from neo4j import GraphDatabase


PENALIZE_HUB = os.environ.get("PENALIZE_HUB", "log").lower()
if PENALIZE_HUB not in {"log", "linear", "none"}:
    raise ValueError(f"PENALIZE_HUB must be one of log|linear|none, got {PENALIZE_HUB!r}")

WEIGHT_EXPR = {
    "log":    "base_weight * (log((toFloat(N_C) + 1.0) / (toFloat(df) + 1.0)) + 1.0)",
    "linear": "CASE WHEN df = 0 THEN base_weight ELSE base_weight / df END",
    "none":   "base_weight",
}[PENALIZE_HUB]


os.makedirs(os.path.expanduser("./matrices/patient_gene/"), exist_ok=True)
driver = GraphDatabase.driver("bolt://localhost:7683", auth=("neo4j", "password"))
with open(os.path.expanduser("./matrices/gene_drug/gene_cols.json")) as f:
    gene_rows = json.load(f)

query = f"""
MATCH (c:Compound) WITH count(c) AS N_C
MATCH (p:Patient)-[r]->(q:Gene)
WITH p, q, N_C,
    CASE r.regulation
        WHEN 'Up'   THEN  1.0
        WHEN 'Down' THEN -1.0
        ELSE              0.0
    END AS base_weight,
    COUNT {{ (:Compound)-[]->(q) }} AS df
RETURN p.name AS patient_name,
       q.name AS gene_name,
       {WEIGHT_EXPR} AS weight
"""

print(f"PENALIZE_HUB={PENALIZE_HUB}")

patient_gene_matrix = []
with driver.session() as session:
    result = session.run(query).to_df().dropna()

for patient_name, df_group in result.groupby("patient_name"):
    gene_name_to_weight = dict(zip(df_group["gene_name"], df_group["weight"]))
    patient_gene_matrix.append(
        [gene_name_to_weight.get(name, 0) for name in gene_rows]
    )

np.save(
    os.path.expanduser("./matrices/patient_gene/patient_gene_matrix.npy"),
    np.array(patient_gene_matrix),
)
