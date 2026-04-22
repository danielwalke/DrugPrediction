import json
import os
with open(os.path.expanduser('./matrices/drug_scores/gene_drugs.json'), 'r') as f:
    gene_drug_scores = json.load(f)

with open(os.path.expanduser('./matrices/drug_scores/protein_drugs.json'), 'r') as f:
    protein_drug_scores = json.load(f)

union_keys = set(gene_drug_scores.keys()) | set(protein_drug_scores.keys())
combined_scores = {}
for key in union_keys:
    gene_score = float(gene_drug_scores.get(key, 0))
    protein_score = float(protein_drug_scores.get(key, 0))
    combined_score = (gene_score + protein_score) / 2
    combined_scores[key] = combined_score


with open(os.path.expanduser('./matrices/drug_scores/combined_drugs.json'), 'w') as f:
    json.dump(dict(sorted(combined_scores.items(), key=lambda item: item[1])), f, indent=4)