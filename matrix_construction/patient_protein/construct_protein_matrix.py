"""
Build the patient x protein regulation matrix.

Per entry:

    R[patient, protein] = sign(regulation) * IDF(protein)

where the per-feature IDF is selected by the environment variable
PENALIZE_HUB (default "log"):

    PENALIZE_HUB=log     (default)  IDF(p) = log( (N_C + 1) / (df_p + 1) ) + 1
    PENALIZE_HUB=linear              IDF(p) = 1 / df_p              (legacy)
    PENALIZE_HUB=none                IDF(p) = 1

Motivation:
    df_p = number of Compounds that target protein p.
    The legacy "linear" weighting 1/df is unstable: a protein hit by
    one compound gets weight 1.0 while a protein hit by 1000 compounds
    gets 0.001, a 1000x dynamic range that rewards obscure targets
    pathologically. The default "log" weighting is the standard
    smoothed IDF from information retrieval and has bounded range
    [1, log N_C + 1] while still monotonically favouring rare targets.
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


os.makedirs(os.path.expanduser("./matrices/patient_protein/"), exist_ok=True)
driver = GraphDatabase.driver("bolt://localhost:7683", auth=("neo4j", "password"))
with open(os.path.expanduser("./matrices/protein_drug/protein_cols.json")) as f:
    protein_rows = json.load(f)

query = f"""
MATCH (c:Compound) WITH count(c) AS N_C
MATCH (p:Patient)-[r]-(q:Protein)
MATCH (p)-[ra]-(g:Gene)--(q) 
WHERE r.regulation = ra.regulation
WITH p, q, N_C,
    CASE r.regulation
        WHEN 'Up'   THEN  1.0
        WHEN 'Down' THEN -1.0
        ELSE              0.0
    END AS base_weight,
    COUNT {{ (:Compound)-[]->(q) }} AS df
RETURN p.name AS patient_name,
       q.name AS protein_name,
       {WEIGHT_EXPR} AS weight
"""

print(f"PENALIZE_HUB={PENALIZE_HUB}")

patient_protein_matrix = []
with driver.session() as session:
    result = session.run(query).to_df().dropna()
print(result.head())
for patient_name, df_group in result.groupby("patient_name"):
    protein_name_to_weight = dict(zip(df_group["protein_name"], df_group["weight"]))
    patient_protein_matrix.append(
        [protein_name_to_weight.get(name, 0) for name in protein_rows]
    )
np.save(
    os.path.expanduser("./matrices/patient_protein/patient_protein_matrix.npy"),
    np.array(patient_protein_matrix),
)
