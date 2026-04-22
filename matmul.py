"""
Patient-drug scoring.

Definition
----------
Let
    R in R^{Q x F}         patient regulation (Up=+1, Down=-1, scaled by the
                           feature-side IDF weight produced by the patient
                           matrix constructor),
    A in {-1, 0, +1}^{F x D}   feature-drug adjacency (protein or gene).
    w in R^{D}             drug-side specificity prior (1 / sqrt(deg(d)+1)).

Score for patient q, drug d:

    S(q, d) = sum_f  R[q, f] * A[f, d] * w[d]
            = (R @ A * w)[q, d].

Interpretation: more negative = better therapeutic match (drug opposes the
patient's dysregulation). This is the convention the rest of the pipeline
(combined_score_analysis.py, debug*.py) already assumes and which
ranking-based optimizers such as scipy.optimize.fmin can minimize directly.

Notes on what this file changed vs. the previous implementation:
- Removed `protein_drug_matrix * sqrt(patient_degree_vector)`: that factor
  is a no-op for a single patient (constant) and dimensionally incoherent
  for multiple patients (it broadcasts a per-patient scalar along the
  feature axis of an F x D matrix). If patient-side normalization is
  desired it belongs on R, not on the drug matrix.
- Added an optional permutation null (`--n-perm`) that complements the
  drug-side hub penalty. The hub penalty is a prior on drug size; the
  permutation null tests whether the *sign pattern* of overlap with the
  patient is non-random (holding the patient's support fixed and
  shuffling only the ± signs). The two are complementary -- see PR
  description for examples.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np


def _resort(items: list, idc: np.ndarray) -> list:
    return [items[int(i)] for i in idc]


def score_all(
    patient_feature_matrix: np.ndarray,   # (Q, F)
    feature_drug_matrix: np.ndarray,      # (F, D)
    drug_weights: np.ndarray,             # (D,)
) -> np.ndarray:
    """Return S with shape (Q, D). More negative = better match."""
    # Apply the drug prior as a column-wise scaling of the F x D matrix,
    # then a single matmul. Broadcasting is aligned on the drug axis.
    weighted = feature_drug_matrix * drug_weights[None, :]    # (F, D)
    return patient_feature_matrix @ weighted                  # (Q, D)


def permutation_z(
    patient_row: np.ndarray,              # (F,)
    feature_drug_matrix: np.ndarray,      # (F, D)
    drug_weights: np.ndarray,             # (D,)
    n_perm: int,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Per-drug z-score against `n_perm` sign permutations of `patient_row`.

    Only the ± signs on the patient's nonzero support are shuffled; the
    support (which features are dysregulated) and the per-feature
    magnitudes are preserved. This isolates sign agreement from raw
    overlap size, which the drug-side hub penalty already handles.

    Returns (observed, z). For this pipeline lower is better, so a very
    *negative* z means the observed score is much more negative than
    expected under random sign shufflings -- strong evidence the drug's
    effect directions align against the patient's dysregulation.
    """
    rng = np.random.default_rng(seed)
    patient_row = patient_row.astype(np.float64, copy=True)
    weighted = feature_drug_matrix * drug_weights[None, :]    # (F, D)

    observed = patient_row @ weighted                         # (D,)

    nz = np.flatnonzero(patient_row)
    magnitudes = np.abs(patient_row[nz])
    signs = np.sign(patient_row[nz]).copy()

    null = np.empty((n_perm, feature_drug_matrix.shape[1]), dtype=np.float64)
    buf = np.zeros_like(patient_row)
    for i in range(n_perm):
        rng.shuffle(signs)
        buf[:] = 0.0
        buf[nz] = signs * magnitudes
        null[i] = buf @ weighted

    mu = null.mean(axis=0)
    sd = null.std(axis=0) + 1e-12
    return observed, (observed - mu) / sd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("patient_feature_matrix")
    p.add_argument("drug_weights")
    p.add_argument("feature_drug_matrix")
    p.add_argument("compound_rows")
    p.add_argument("output_drug_scores")
    p.add_argument(
        "--n-perm",
        type=int,
        default=0,
        help="If >0, also compute a permutation-null z-score per drug "
        "and write it to <output>.z.json. Default: 0 (off).",
    )
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(os.path.expanduser("./matrices/drug_scores"), exist_ok=True)

    patient_feature_matrix = np.load(os.path.expanduser(args.patient_feature_matrix))
    drug_weights = np.load(os.path.expanduser(args.drug_weights))
    feature_drug_matrix = np.load(os.path.expanduser(args.feature_drug_matrix)).T

    print("Patient-feature matrix:", patient_feature_matrix.shape)
    print("Feature-drug matrix:   ", feature_drug_matrix.shape)
    print("Drug weights:          ", drug_weights.shape)

    S = score_all(patient_feature_matrix, feature_drug_matrix, drug_weights)
    print("Patient-drug score matrix:", S.shape)

    # One patient for now; this indexes row 0 deliberately. See README.
    patient_idx = 0
    scores = S[patient_idx]
    sorted_idx = np.argsort(scores)   # ascending -> most negative first

    with open(os.path.expanduser(args.compound_rows)) as f:
        compound_rows = json.load(f)

    drug_scores = dict(
        zip(_resort(compound_rows, sorted_idx), scores[sorted_idx].astype(float).tolist())
    )
    with open(os.path.expanduser(args.output_drug_scores), "w") as f:
        json.dump(drug_scores, f, indent=4)

    for i, idx in enumerate(sorted_idx):
        if i == 10:
            break
        print(
            f"  #{i:<3d} {compound_rows[idx]:<60s} "
            f"score={scores[idx]:+.4f}"
        )

    if args.n_perm > 0:
        observed, z = permutation_z(
            patient_feature_matrix[patient_idx],
            feature_drug_matrix,
            drug_weights,
            n_perm=args.n_perm,
            seed=args.seed,
        )
        z_sorted = np.argsort(z)      # ascending: most negative = strongest
        z_path = args.output_drug_scores.rstrip(".json") + ".z.json"
        z_out = {
            compound_rows[int(i)]: {
                "score": float(observed[i]),
                "zscore": float(z[i]),
            }
            for i in z_sorted
        }
        with open(os.path.expanduser(z_path), "w") as f:
            json.dump(z_out, f, indent=4)
        print(f"\nWrote permutation z-scores ({args.n_perm} perms) to {z_path}")
        print("Top 10 by |z| (most negative):")
        for i, idx in enumerate(z_sorted[:10]):
            print(
                f"  #{i:<3d} {compound_rows[int(idx)]:<60s} "
                f"score={observed[idx]:+.4f}  z={z[idx]:+.3f}"
            )


if __name__ == "__main__":
    main()
