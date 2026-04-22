"""
Combine protein-based and gene-based drug scores.

Why we need something better than a simple mean
-----------------------------------------------
The previous implementation averaged the two raw scores:

    S_combined(d) = (S_gene(d) + S_protein(d)) / 2.

The problem: the two channels live on different scales.

* The protein graph has ~14 relation types and many more edges per drug
  than the gene graph (which has 3 relation types). The spread
  (standard deviation) of S_protein is therefore typically much larger
  than S_gene.
* When you average two quantities with different spreads, the one with
  the bigger spread dominates. A drug that is mediocre in proteins but
  excellent in genes will be outranked by a drug that is merely
  average in proteins, because "average in proteins" already has more
  absolute magnitude than "excellent in genes".

We provide two principled combiners. Both treat the two channels
symmetrically.

1. Z-score average
------------------
    tilde S_c(d) = ( S_c(d) - mean_c ) / std_c     for c in {gene, protein}
    S_zcomb(d)   = 0.5 * ( tilde S_gene(d) + tilde S_protein(d) )

Read "tilde S" as "how many standard deviations away from the
channel's mean". Averaging these puts both channels on the same
footing. Lower S_zcomb = better (same sign convention as the raw
scores in this pipeline).

2. Reciprocal Rank Fusion (RRF)
-------------------------------
    rank_c(d) = 1 for the drug with the most negative S_c (best),
                2 for the second best, and so on.
    S_rrf(d) = sum_c  1 / (k + rank_c(d)),           k = 60.

This is the standard combiner in information retrieval. It only uses
the *order* of drugs in each channel -- it is invariant to the raw
scale of S_gene and S_protein, and it is insensitive to outliers. The
constant k = 60 is the IR default and controls how much weight is
given to drugs that are top-ranked in just one channel. Higher
S_rrf = better.

RRF is our recommended ranking criterion; the z-score is kept because
it gives a physical number ("2 sigma better than average").
"""
from __future__ import annotations

import json
import os
from typing import Dict

import numpy as np


def zscore(values: Dict[str, float]) -> Dict[str, float]:
    v = np.asarray(list(values.values()), dtype=np.float64)
    mu, sd = v.mean(), v.std() + 1e-12
    return {k: (x - mu) / sd for k, x in values.items()}


def rrf_ranks(values: Dict[str, float], k: int = 60) -> Dict[str, float]:
    """Most-negative raw score gets rank 1 (best)."""
    ordered = sorted(values.items(), key=lambda kv: kv[1])
    return {key: 1.0 / (k + rank) for rank, (key, _) in enumerate(ordered, start=1)}


def main(alpha: float = 0.5) -> None:
    with open(os.path.expanduser("./matrices/drug_scores/gene_drugs.json")) as f:
        gene = {k: float(v) for k, v in json.load(f).items()}
    with open(os.path.expanduser("./matrices/drug_scores/protein_drugs.json")) as f:
        prot = {k: float(v) for k, v in json.load(f).items()}

    gene_z, prot_z = zscore(gene), zscore(prot)
    gene_r, prot_r = rrf_ranks(gene), rrf_ranks(prot)

    union = set(gene) | set(prot)
    combined: Dict[str, Dict[str, float]] = {}
    for key in union:
        gz = gene_z.get(key, 0.0)
        pz = prot_z.get(key, 0.0)
        gr = gene_r.get(key, 0.0)
        pr = prot_r.get(key, 0.0)
        combined[key] = {
            # Lower = better, same convention as the raw scores.
            "zscore_combined": alpha * pz + (1.0 - alpha) * gz,
            # Higher = better (this is the recommended ranking).
            "rrf_combined": gr + pr,
            "protein_score": prot.get(key, 0.0),
            "gene_score": gene.get(key, 0.0),
        }

    # Write two outputs for convenience:
    # - combined_drugs.json: sorted ascending by zscore_combined (lowest = best),
    #   preserving the legacy "lowest = best" convention.
    # - combined_drugs_rrf.json: sorted descending by rrf_combined (highest = best).
    out_dir = os.path.expanduser("./matrices/drug_scores")
    os.makedirs(out_dir, exist_ok=True)

    by_z = dict(sorted(combined.items(), key=lambda kv: kv[1]["zscore_combined"]))
    with open(os.path.join(out_dir, "combined_drugs.json"), "w") as f:
        json.dump(by_z, f, indent=4)

    by_rrf = dict(sorted(combined.items(), key=lambda kv: -kv[1]["rrf_combined"]))
    with open(os.path.join(out_dir, "combined_drugs_rrf.json"), "w") as f:
        json.dump(by_rrf, f, indent=4)

    print("Top 10 by RRF (recommended):")
    for i, (name, row) in enumerate(by_rrf.items()):
        if i >= 10:
            break
        print(
            f"  {name:<60s}  rrf={row['rrf_combined']:.4f}  "
            f"z={row['zscore_combined']:+.3f}"
        )


if __name__ == "__main__":
    main()
