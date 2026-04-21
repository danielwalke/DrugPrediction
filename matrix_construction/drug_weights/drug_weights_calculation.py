import json
import numpy as np
import pandas as pd
import os

os.makedirs(os.path.expanduser('./matrices/drug_weights/'), exist_ok=True)

with open(os.path.expanduser('./matrices/drug_weights/compound_metrics.json'), 'r') as f:
    compound_metrics = json.load(f)

with open(os.path.expanduser('./matrices/protein_drug/compound_rows.json'), 'r') as f:
    compound_rows = json.load(f)

weights = []
for comp in compound_rows:
    weight = 1 / (1 + compound_metrics[comp]["out_degree"])
    weights.append(weight)
weights = np.array(weights)
weights = weights / np.sum(weights)
print(pd.DataFrame({'compound': compound_rows, 'weight': weights}).describe())
np.save(os.path.expanduser('./matrices/drug_weights/drug_weights.npy'), weights)

