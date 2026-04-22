#!/usr/bin/env bash
# Run the full drug-prediction pipeline.
#
# Optional knobs (set as env vars before calling, e.g.
#   PENALIZE_HUB=log DEGREE_SCOPE=targets N_PERM=200 ./drug_prediction.sh):
#
#   PENALIZE_HUB = log (default) | linear | none
#       Per-feature IDF weighting on the patient vectors.
#         log    :  w_f = log((N_C + 1) / (df_f + 1)) + 1   (recommended)
#         linear :  w_f = 1 / df_f                          (legacy, unstable)
#         none   :  w_f = 1
#
#   DEGREE_SCOPE = targets (default) | all
#       Which edges contribute to the drug-side specificity prior
#       w_d = 1 / sqrt(out_degree(d) + 1).
#         targets :  only edges to :Compound, :Protein, :Gene
#         all     :  every label (legacy apoc.node.degree)
#
#   N_PERM = 0 (default) | positive integer
#       If >0, matmul.py also writes a permutation-null z-score per
#       drug to <output>.z.json, shuffling only the ± signs on the
#       patient's support. Useful once pipelines are set up.
#
set -euo pipefail

: "${PENALIZE_HUB:=log}"
: "${DEGREE_SCOPE:=targets}"
: "${N_PERM:=0}"

NPERM_ARGS=()
if [[ "$N_PERM" -gt 0 ]]; then
    NPERM_ARGS=(--n-perm "$N_PERM")
fi

export PENALIZE_HUB DEGREE_SCOPE

python matrix_construction/protein_drug/create_drug_effect_df.py
python matrix_construction/drug_weights/degree-calculation.py \
    "./matrices/protein_drug/compound_rows.json" \
    "./matrices/drug_weights/compound_metrics.json"
python matrix_construction/drug_weights/drug_weights_calculation.py \
    "./matrices/drug_weights/compound_metrics.json" \
    "./matrices/protein_drug/compound_rows.json" \
    "./matrices/drug_weights/drug_weights.npy"
python matrix_construction/patient_protein/construct_protein_matrix.py
python matmul.py \
    "./matrices/patient_protein/patient_protein_matrix.npy" \
    "./matrices/drug_weights/drug_weights.npy" \
    "./matrices/protein_drug/dense_adjacency_matrix.npy" \
    "./matrices/protein_drug/compound_rows.json" \
    "./matrices/drug_scores/protein_drugs.json" \
    "${NPERM_ARGS[@]}"

python matrix_construction/gene_drug/create_drug_effect_df.py
python matrix_construction/drug_weights/degree-calculation.py \
    "./matrices/gene_drug/compound_rows.json" \
    "./matrices/drug_weights/compound_metrics_genes.json"
python matrix_construction/drug_weights/drug_weights_calculation.py \
    "./matrices/drug_weights/compound_metrics_genes.json" \
    "./matrices/gene_drug/compound_rows.json" \
    "./matrices/drug_weights/drug_weights_genes.npy"
python matrix_construction/gene_patient/construct_gene_matrix.py
python matmul.py \
    "./matrices/patient_gene/patient_gene_matrix.npy" \
    "./matrices/drug_weights/drug_weights_genes.npy" \
    "./matrices/gene_drug/dense_adjacency_matrix.npy" \
    "./matrices/gene_drug/compound_rows.json" \
    "./matrices/drug_scores/gene_drugs.json" \
    "${NPERM_ARGS[@]}"

python combined_score_analysis.py
