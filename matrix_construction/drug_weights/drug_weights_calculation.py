import json
import numpy as np
import pandas as pd
import os
import sys

if len(sys.argv) != 4:
    print("Usage: python drug_weights_calculation.py <input_file_path_compound_metrics> <input_file_path_compound_rows> <output_file_path>")
    sys.exit(1)

input_file_path_compound_metrics = sys.argv[1]
input_file_path_compound_rows = sys.argv[2]
output_file_path = sys.argv[3]

os.makedirs(os.path.expanduser('./matrices/drug_weights/'), exist_ok=True)

with open(os.path.expanduser(input_file_path_compound_metrics), 'r') as f:
    compound_metrics = json.load(f)

with open(os.path.expanduser(input_file_path_compound_rows), 'r') as f:
    compound_rows = json.load(f)

weights = []
for comp in compound_rows:
    weight = 1 / np.sqrt(compound_metrics[comp]["out_degree"] + 1)
    weights.append(weight)
weights = np.array(weights)
weights = weights
print(pd.DataFrame({'compound': compound_rows, 'weight': weights}).describe())
np.save(os.path.expanduser(output_file_path), weights)

