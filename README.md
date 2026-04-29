# DrugPrediction

DrugPrediction is a Python-based repository designed to predict drug efficacy for patients using both protein-based and gene-based prediction models. It leverages graph and network data from [PharMeBINet](https://pharmebi.net/) to construct patient, protein, gene, and drug matrices. By applying matrix multiplication algorithms along with degree-based drug weighting, the tool mathematically determines the most effective drugs for a given profile and combines the scores for comprehensive analysis.

## Table of Contents

* [Prerequisites](#prerequisites)
* [Setup and Installation](#setup-and-installation)
* [Protein-Based Prediction](#protein-based-prediction)
* [Gene-Based Prediction](#gene-based-prediction)
* [Combined Analysis](#combined-analysis)

## Prerequisites

Before setting up the project, ensure you have downloaded the latest PharMeBINet database dump:

* [PharMeBINet Download](https://pharmebi.net/#/download)

## Setup and Installation

1. Create the database dump:

```bash
sh pharmebinet/create_db.sh
```

2. Start the database:

```bash
sh pharmebinet/start_db.sh
```

3. Set up the Python virtual environment and install the required dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Protein-Based Prediction

Execute the following commands in sequence to run the protein-based prediction pipeline:

1. Construct the protein-drug matrix (PxD matrix):

```bash
python matrix_construction/protein_drug/create_drug_effect_df.py
```

2. Create the drug degree JSON:

```bash
python matrix_construction/drug_weights/degree-calculation.py "./matrices/protein_drug/compound_rows.json" "./matrices/drug_weights/compound_metrics.json"
```

3. Calculate the drug weights based on these degrees:

```bash
python matrix_construction/drug_weights/drug_weights_calculation.py "./matrices/drug_weights/compound_metrics.json" "./matrices/protein_drug/compound_rows.json" "./matrices/drug_weights/drug_weights.npy"
```

4. Construct the protein-to-patient matrix:

```bash
python matrix_construction/patient_protein/construct_protein_matrix.py
```

5. Use matrix multiplication to extract the highest-scoring drugs:

```bash
python matmul.py "./matrices/patient_protein/patient_protein_matrix.npy" "./matrices/drug_weights/drug_weights.npy" "./matrices/protein_drug/dense_adjacency_matrix.npy" "./matrices/protein_drug/compound_rows.json" "./matrices/drug_scores/protein_drugs.json"
```

## Gene-Based Prediction

Execute the following commands in sequence to run the gene-based prediction pipeline:

1. Construct the gene-drug matrix (PxD matrix):

```bash
python matrix_construction/gene_drug/create_drug_effect_df.py
```

2. Create the drug degree matrix:

```bash
python matrix_construction/drug_weights/degree-calculation.py "./matrices/gene_drug/compound_rows.json" "./matrices/drug_weights/compound_metrics_genes.json"
```

3. Calculate the drug weights based on these degrees:

```bash
python matrix_construction/drug_weights/drug_weights_calculation.py "./matrices/drug_weights/compound_metrics_genes.json" "./matrices/gene_drug/compound_rows.json" "./matrices/drug_weights/drug_weights_genes.npy"
```

4. Construct the gene-to-patient matrix:

```bash
python matrix_construction/gene_patient/construct_gene_matrix.py
```

5. Use matrix multiplication to extract the highest-scoring drugs:

```bash
python matmul.py "./matrices/patient_gene/patient_gene_matrix.npy" "./matrices/drug_weights/drug_weights_genes.npy" "./matrices/gene_drug/dense_adjacency_matrix.npy" "./matrices/gene_drug/compound_rows.json" "./matrices/drug_scores/gene_drugs.json"
```

## Combined Analysis

Once both the protein-based and gene-based prediction pipelines have been executed, you can combine the scores for a final, consolidated analysis.

Run the combination script:

```bash
python combined_score_analysis.py
```

This writes two files to `./matrices/drug_scores/`:

* `combined_drugs.json` — sorted ascending by `zscore_combined` (lowest = best, same sign convention as the raw channel scores).
* `combined_drugs_rrf.json` — sorted descending by `rrf_combined` (highest = best). Recommended as the final ranking.

Each drug entry looks like:

```json
"SomeDrugName": {
    "zscore_combined": -4.5578,
    "rrf_combined":     0.01645,
    "protein_score":   -3.1021,
    "gene_score":      -0.8842
}
```

The four fields, in plain terms:

* `protein_score` — the raw score `S(q, d)` from the protein pipeline (sum of signed, IDF-weighted feature agreement, scaled by the drug specificity prior). Lower = better.
* `gene_score` — the raw score `S(q, d)` from the gene pipeline. Lower = better.
* `zscore_combined` — the two raw scores put on a common scale, then averaged. Lower = better.
* `rrf_combined` — the two raw scores converted to ranks, then fused. Higher = better.

### Why two combined scores?

The two raw channels live on different scales. The protein graph has ~14 distinct relation types and is much denser than the gene graph (~3 relation types), so `std(protein_score) ≫ std(gene_score)` in practice. If you just averaged them, the protein channel would dominate: a drug that is average in proteins but excellent in genes would lose to a drug that is merely average in proteins, because "average in proteins" already has more absolute magnitude than "excellent in genes".

`zscore_combined` and `rrf_combined` are two principled ways to put the channels on equal footing. They almost always agree at the top of the ranking, but they weight the tail differently.

### `zscore_combined` — standardise, then average

For each channel `c ∈ {protein, gene}`, compute the mean `μ_c` and standard deviation `σ_c` of that channel's raw scores across all drugs. Standardise:

```
z_c(d) = (S_c(d) − μ_c) / σ_c
```

`z_c(d)` is the number of standard deviations drug `d` is below (or above) the channel mean. Because each channel has been standardised, both channels contribute the same spread to the combined score. Then average:

```
zscore_combined(d) = 0.5 · z_protein(d) + 0.5 · z_gene(d)
```

Interpretation: "this drug is on average this many σ below the mean of each channel". `-4.56` means roughly 4.56 standard deviations below the mean when averaged across the two channels — strong, interpretable evidence that `d` is unusually good at opposing the patient's dysregulation in both modalities.

Sign convention: lower = better, matching the raw scores.

### `rrf_combined` — Reciprocal Rank Fusion

For each channel `c`, rank the drugs from best (rank 1, the most-negative raw score) to worst. Then:

```
rrf_combined(d) = Σ_c  1 / (k + rank_c(d)),        k = 60
```

The constant `k = 60` is the standard choice in information retrieval. It controls how sharply the top of each ranking is rewarded: with `k = 60`, rank 1 contributes `1/61 ≈ 0.0164`, rank 2 contributes `1/62 ≈ 0.0161`, rank 100 contributes `1/160 ≈ 0.0063`, rank 10000 contributes `1/10060 ≈ 0.0001`. The maximum achievable RRF is `2 / 61 ≈ 0.0328` (rank 1 in both channels). `0.01645` in the example above means the drug is rank 1 in one channel and very low in the other, or middle-of-the-pack in both.

RRF uses only the *order* of drugs, not the raw magnitudes. That makes it:

* **Scale-free** — irrelevant how wide each channel's distribution is.
* **Outlier-robust** — one channel having a single pathological score doesn't shift the ranking.
* **Less informative** when you want to talk about "how much better than average" a drug is (use `zscore_combined` for that).

Sign convention: higher = better.

### Which should I use?

The recommended ranking is `rrf_combined` — it is the standard IR combiner and is robust to the scale differences between the protein and gene channels. `zscore_combined` is kept because it carries a physical interpretation ("σ below the mean") that is easier to communicate.

In practice the top-10s of the two files will overlap heavily; divergences happen further down the list and are a useful signal for case-by-case review.

See the docstring at the top of `combined_score_analysis.py` for the exact implementation.

## Configuration knobs

The pipeline exposes three optional environment variables. All have sensible defaults, so `sh drug_prediction.sh` works unchanged.

| Variable | Default | Values | Effect |
| --- | --- | --- | --- |
| `PENALIZE_HUB` | `log` | `log`, `linear`, `none` | Per-feature IDF weighting on the patient vectors. `log` is the standard smoothed IDF `log((N_C+1)/(df+1))+1`. `linear` reproduces the legacy `1/df`. `none` disables it. |
| `DEGREE_SCOPE` | `targets` | `targets`, `all` | Which edges count towards the drug specificity prior `w_d = 1/√(out_degree+1)`. `targets` counts only edges to `:Compound`, `:Protein` and `:Gene` (recommended). `all` counts every label. |
| `N_PERM` | `0` | non-negative integer | If > 0, `matmul.py` also writes a permutation-null z-score per drug to `<output>.z.json`. Shuffles only the ± signs on the patient's support, holding the support fixed. Complements the drug-size prior — large negative z means "sign alignment with this drug is much better than chance". 200–1000 is typical. |

Example:

```bash
PENALIZE_HUB=log DEGREE_SCOPE=targets N_PERM=500 sh drug_prediction.sh
```

Investigate Neo4j Results
MATCH path=(d:Compound)-[]-(p)-[]-(n:Patient) 
WHERE d.name = "Dabrafenib" AND ("Gene" IN labels(p) OR "Protein" IN labels(p)) 
RETURN path


MATCH path=(d:Compound)-[]-(p)-[]-(n:Patient) 
WHERE d.name = "Pyrophosphate 2-" AND ("Gene" IN labels(p) OR "Protein" IN labels(p)) 
RETURN path


## Scoring params
1) Default not: Side effect penalization with patient_normal_feature_matrix = (patient_feature_matrix == 0).astype(int) * 0 
2) Default yes: Limiting  number of affected genes and proteins to only those where gene and proteins show similar expression trends: "MATCH (p:Patient)-[r]->(q:Gene) MATCH (p)-[ra]-(a:Protein)--(q) WHERE r.regulation = ra.regulation"
3) Default yes: Down-weighting drugs with many interactions to other compounds and chemicals besides RESEMBLES interactions
4) Default not: Weighting the influence of dysregulated proteins and genes higher than the physiolgoically expressed ones