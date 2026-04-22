#!/usr/bin/env bash
# Run the full drug-prediction pipeline.
#
# =============================================================================
# Symbol legend (same symbols used across matmul.py, the matrix constructors
# and this script):
#
#   q         a patient (one row of the patient matrix).
#   d         a compound (a candidate drug, one column of the drug matrix).
#   f         a feature: a protein in the protein pipeline, or a gene in the
#             gene pipeline.
#   N_C       total number of :Compound nodes in PharMeBINet. Used as the
#             "corpus size" in the IDF formula below.
#   df_f      "document frequency" of feature f, i.e. the number of
#             compounds that target f (number of edges :Compound -[]-> f).
#   R[q, f]   signed regulation of feature f for patient q, in {-1, 0, +1},
#             multiplied by the per-feature IDF weight w_f.
#   A[f, d]   effect of drug d on feature f, in {-1, 0, +1}. Upregulation,
#             agonism = +1; inhibition, degradation = -1; unknown = 0
#             (see config/compound_{,gene_}connectivity.json).
#   out_degree(d)  number of outgoing edges from compound d (scope below).
#   w_d       drug-side specificity prior, 1 / sqrt(out_degree(d) + 1).
#             Penalises compounds that interact with many entities.
#   S(q, d)   patient-drug score = sum_f R[q, f] * A[f, d] * w_d.
#             More negative = better therapeutic match (drug opposes the
#             patient's dysregulation).
# =============================================================================
#
# Optional knobs (set as env vars before calling, e.g.
#   PENALIZE_HUB=log DEGREE_SCOPE=targets N_PERM=500 ./drug_prediction.sh):
#
# -----------------------------------------------------------------------------
# PENALIZE_HUB = log (default) | linear | none
# -----------------------------------------------------------------------------
#   Controls w_f, the per-feature IDF weight applied inside R[q, f]. The
#   IDF weight down-weights features that are hit by many compounds (because
#   they carry less information about which specific drug fits the patient),
#   and up-weights features that are hit by few.
#
#     log    :  w_f = log((N_C + 1) / (df_f + 1)) + 1      (recommended)
#               Smoothed inverse-document-frequency, the standard choice
#               from information retrieval. Bounded in [1, log N_C + 1].
#
#     linear :  w_f = 1 / df_f                             (legacy)
#               The original formula. A feature hit by 1 compound gets
#               weight 1.0; a feature hit by 1000 compounds gets 0.001.
#               The 1000x dynamic range is unstable and over-rewards
#               obscure targets.
#
#     none   :  w_f = 1
#               No IDF weighting; treat all features equally.
#
# -----------------------------------------------------------------------------
# DEGREE_SCOPE = targets (default) | all
# -----------------------------------------------------------------------------
#   Controls which edges count towards out_degree(d) in the drug-side
#   specificity prior w_d = 1 / sqrt(out_degree(d) + 1).
#
#     targets : only edges to :Compound, :Protein and :Gene. This is
#               the "interaction promiscuity" scope: a drug that binds
#               many proteins, regulates many genes or has many drug-drug
#               interactions is penalised. Edges to :Disease, :SideEffect,
#               :Pathway, :Anatomy etc. are excluded because a drug that
#               "treats many diseases" is not therefore unspecific.
#
#     all     : every outgoing edge (legacy apoc.node.degree(c)).
#
# -----------------------------------------------------------------------------
# N_PERM = 0 (default) | positive integer (200-1000 typical)
# -----------------------------------------------------------------------------
#   If >0, matmul.py also writes a permutation-null z-score per drug to
#   <output>.z.json. Procedure:
#
#     1. Fix the patient's support (which features are nonzero in R[q, :]).
#        Record the magnitudes |R[q, f]| on the support.
#     2. Take the observed score S(q, d) for each drug d.
#     3. Repeat N_PERM times: shuffle only the +/- signs on the support,
#        leaving the support and magnitudes unchanged, and recompute the
#        score against every drug.
#     4. For each drug d, fit a mean mu_d and std sigma_d to the N_PERM
#        shuffled scores and report
#               z(q, d) = (S(q, d) - mu_d) / sigma_d.
#
#   Large negative z means the observed alignment of signs between the
#   patient's dysregulation and drug d is far more anti-correlated than
#   random sign assignments would produce. That is evidence that the
#   direction of d's effect genuinely opposes the patient's state --
#   not merely that d happens to touch many of the same features.
#
#   Worked example:
#     Patient has 20 dysregulated proteins (R[q, :] nonzero at 20 positions).
#     Drug D hits 500 proteins total, of which all 20 of the patient's
#     dysregulated proteins are included and every one of D's effects on
#     those 20 is the opposite sign of the patient's regulation.
#     * The drug-size prior w_D = 1/sqrt(501) heavily damps D.
#     * The permutation null reshuffles the 20 +/- signs of the patient.
#       Only the original sign pattern produces the observed very-negative
#       S(q, D); random sign flips on the same 20 positions almost never
#       match D's effect directions simultaneously on all 20.
#       -> mu_D is close to zero, sigma_D small, z_D is a large negative
#          number. D is promoted in the zscore ranking even though it
#          was damped in the raw ranking.
#
#   Output files are written next to the raw scores:
#     ./matrices/drug_scores/protein_drugs.z.json
#     ./matrices/drug_scores/gene_drugs.z.json
#   Each entry is {"score": S(q, d), "zscore": z(q, d)}. The JSON is
#   sorted by most-negative z first.
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
