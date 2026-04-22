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

* `combined_drugs.json` — sorted ascending by z-score-averaged score (lowest = best, same convention as the raw channel scores).
* `combined_drugs_rrf.json` — sorted descending by Reciprocal Rank Fusion (highest = best). Recommended as the final ranking.

See the docstring at the top of `combined_score_analysis.py` for the math.

## Configuration knobs

The pipeline exposes three optional environment variables. All have sensible defaults, so `sh drug_prediction.sh` works unchanged.

| Variable | Default | Values | Effect |
| --- | --- | --- | --- |
| `PENALIZE_HUB` | `log` | `log`, `linear`, `none` | Per-feature IDF weighting on the patient vectors. `log` is the standard smoothed IDF `log((N_C+1)/(df+1))+1`. `linear` reproduces the legacy `1/df`. `none` disables it. |
| `DEGREE_SCOPE` | `targets` | `targets`, `all` | Which edges count towards the drug specificity prior `w_d = 1/√(out_degree+1)`. `targets` counts only edges to `:Compound`, `:Protein` and `:Gene` (recommended). `all` counts every label (legacy). |
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