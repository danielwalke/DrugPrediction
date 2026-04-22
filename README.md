# DrugPrediction

DrugPrediction is a Python-based repository designed to predict drug efficacy for patients using both protein-based and gene-based prediction models. It leverages graph and network data from [PharMeBINet](https://www.google.com/search?q=https://pharmebi.net/) to construct patient, protein, gene, and drug matrices. By applying matrix multiplication algorithms along with degree-based drug weighting, the tool mathematically determines the most effective drugs for a given profile and combines the scores for comprehensive analysis.

## Table of Contents

* [Prerequisites](https://www.google.com/search?q=%23prerequisites)

* [Setup and Installation](https://www.google.com/search?q=%23setup-and-installation)

* [Protein-Based Prediction](https://www.google.com/search?q=%23protein-based-prediction)

* [Gene-Based Prediction](https://www.google.com/search?q=%23gene-based-prediction)

* [Combined Analysis](https://www.google.com/search?q=%23combined-analysis)

## Prerequisites

Before setting up the project, ensure you have downloaded the latest PharMeBINet database dump:

* [PharMeBINet Download](https://pharmebi.net/#/download)

## Setup and Installation

1. Create the database dump:

```
sh pharmebinet/create_db.sh

```

2. Start the database:

```
sh pharmebinet/start_db.sh

```

3. Set up the Python virtual environment and install the required dependencies:

```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

```

## Protein-Based Prediction

Execute the following commands in sequence to run the protein-based prediction pipeline:

1. Construct the protein-drug matrix (PxD matrix):

```
python matrix_construction/protein_drug/create_drug_effect_df.py

```

2. Create the drug degree JSON:

```
python matrix_construction/drug_weights/degree-calculation.py "./matrices/protein_drug/compound_rows.json" "./matrices/drug_weights/compound_metrics.json"

```

3. Calculate the drug weights based on these degrees:

```
python matrix_construction/drug_weights/drug_weights_calculation.py "./matrices/drug_weights/compound_metrics.json" "./matrices/protein_drug/compound_rows.json" "./matrices/drug_weights/drug_weights.npy"

```

4. Construct the protein-to-patient matrix:

```
python matrix_construction/patient_protein/construct_protein_matrix.py

```

5. Use matrix multiplication to extract the highest-scoring drugs:

```
python matmul.py "./matrices/patient_protein/patient_protein_matrix.npy" "./matrices/drug_weights/drug_weights.npy" "./matrices/protein_drug/dense_adjacency_matrix.npy" "./matrices/protein_drug/compound_rows.json" "./matrices/drug_scores/protein_drugs.json"

```

## Gene-Based Prediction

Execute the following commands in sequence to run the gene-based prediction pipeline:

1. Construct the gene-drug matrix (PxD matrix):

```
python matrix_construction/gene_drug/create_drug_effect_df.py

```

2. Create the drug degree matrix:

```
python matrix_construction/drug_weights/degree-calculation.py "./matrices/gene_drug/compound_rows.json" "./matrices/drug_weights/compound_metrics_genes.json"

```

3. Calculate the drug weights based on these degrees:

```
python matrix_construction/drug_weights/drug_weights_calculation.py "./matrices/drug_weights/compound_metrics_genes.json" "./matrices/gene_drug/compound_rows.json" "./matrices/drug_weights/drug_weights_genes.npy"

```

4. Construct the gene-to-patient matrix:

```
python matrix_construction/gene_patient/construct_gene_matrix.py

```

5. Use matrix multiplication to extract the highest-scoring drugs:

```
python matmul.py "./matrices/patient_gene/patient_gene_matrix.npy" "./matrices/drug_weights/drug_weights_genes.npy" "./matrices/gene_drug/dense_adjacency_matrix.npy" "./matrices/gene_drug/compound_rows.json" "./matrices/drug_scores/gene_drugs.json"

```

## Combined Analysis

Once both the protein-based and gene-based prediction pipelines have been executed, you can combine the scores for a final, consolidated analysis.

Run the combination script:

```
python combined_score_analysis.py

```

This will output a JSON file containing the combined results, located at:
`./matrices/drug_scores/combined_drugs.json`