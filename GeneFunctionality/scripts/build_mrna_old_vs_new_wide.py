#!/usr/bin/env python3
"""
Build the wide-format CSV that get_mrna_transcript_coords.py and
extract_gwas_data_mrna.sh expect, using the mrna_old_vs_new feature CSVs and
mRNA model predictions from the predicting-genic-features sibling project.

Pairs exon2 and exon3 by row index (verified: row N of v7-exon2 == row N of
v7-exon3 by GeneID, same for v49). Concatenates v7 (Functional=No) and v49
(Functional=Yes) into a single 500-row wide CSV.

Output column layout matches data/model_predictions/matrix-positive-100mrna_
exon2-exon3-features.csv (39 cols, no transcript columns yet — those are
added by get_mrna_transcript_coords.py downstream).
"""

import argparse
import os
import sys
import pandas as pd

PGF = "../predicting-genic-features"
DEFAULT_FEATURES_DIR = f"{PGF}/results/mrna_old_vs_new"
DEFAULT_DATASETS_DIR = f"{PGF}/data/datasets/mrna_old_vs_new"

SAMPLES = ["v7-exon2", "v7-exon3", "v49-exon2", "v49-exon3"]
FEATURE_COLS = [
    "GC_percentage", "phyloP_max_241w", "phyloP_max_100w",
    "RPKM_tissue", "RPKM_primary-cell", "copy_number", "repeat_distance",
    "Interaction_ave", "coding_potential", "Max_covariance", "MFE", "Random",
]


def build_pair(source: str, features_dir: str, datasets_dir: str) -> pd.DataFrame:
    """Pair source-exon2 with source-exon3 by row index.

    Real gene symbols come from the source dataset CSVs (the features CSVs
    have their GeneID column overwritten with ID by merge_features.py).
    The mRNA prediction CSV is the combined 1000-row file (concat order:
    v7-exon2, v7-exon3, v49-exon2, v49-exon3 — see
    scripts/run_mrna_old_vs_new_merge_and_predict.sh Step 4).
    """
    ex2_feats = pd.read_csv(f"{features_dir}/{source}-exon2_features.csv").reset_index(drop=True)
    ex3_feats = pd.read_csv(f"{features_dir}/{source}-exon3_features.csv").reset_index(drop=True)
    ex2_src = pd.read_csv(f"{datasets_dir}/{source}-exon2-dataset.csv").reset_index(drop=True)
    ex3_src = pd.read_csv(f"{datasets_dir}/{source}-exon3-dataset.csv").reset_index(drop=True)
    all_preds = pd.read_csv(f"{features_dir}/predictions/mrna_model_predictions.csv").reset_index(drop=True)

    offsets = {"v7-exon2": 0, "v7-exon3": 250, "v49-exon2": 500, "v49-exon3": 750}
    ex2_preds = all_preds.iloc[offsets[f"{source}-exon2"]:offsets[f"{source}-exon2"]+250].reset_index(drop=True)
    ex3_preds = all_preds.iloc[offsets[f"{source}-exon3"]:offsets[f"{source}-exon3"]+250].reset_index(drop=True)

    # Sanity: gene-symbol alignment between exon2 and exon3 in the source CSVs
    if not (ex2_src["GeneID"].values == ex3_src["GeneID"].values).all():
        sys.exit(f"GeneID mismatch for {source} ex2/ex3 pairing")

    rows = []
    for i in range(len(ex2_feats)):
        e2 = ex2_feats.iloc[i]
        e3 = ex3_feats.iloc[i]
        s2 = ex2_src.iloc[i]
        p2 = ex2_preds.iloc[i]
        p3 = ex3_preds.iloc[i]

        row = {
            "ID": f"{source}_gene_{i+1}",
            "ex2_Chromosome": e2["Chromosome"],
            "ex2_Start": int(e2["Start"]),
            "ex2_End": int(e2["End"]),
            "ex2_Sequence": e2["Sequence"],
            "ex2_prob_No": p2["prob_No"],
            "ex2_prob_Yes": p2["prob_Yes"],
            "ex3_Chromosome": e3["Chromosome"],
            "ex3_Start": int(e3["Start"]),
            "ex3_End": int(e3["End"]),
            "ex3_Sequence": e3["Sequence"],
            "ex3_prob_No": p3["prob_No"],
            "ex3_prob_Yes": p3["prob_Yes"],
            "GeneID": s2["GeneID"],
            # Functional uses Yes/No to match existing wide-format files.
            "Functional": "Yes" if int(e2["Functional"]) == 1 else "No",
        }
        for fc in FEATURE_COLS:
            row[f"ex2_{fc}"] = e2[fc]
            row[f"ex3_{fc}"] = e3[fc]
        rows.append(row)

    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features-dir", default=DEFAULT_FEATURES_DIR,
                    help="Directory containing the *_features.csv and predictions/")
    ap.add_argument("--datasets-dir", default=DEFAULT_DATASETS_DIR,
                    help="Directory containing the original *-dataset.csv files (for real gene symbols)")
    ap.add_argument("--output", required=True,
                    help="Output wide CSV (no transcript columns yet)")
    args = ap.parse_args()

    v7 = build_pair("v7", args.features_dir, args.datasets_dir)
    v49 = build_pair("v49", args.features_dir, args.datasets_dir)
    print(f"v7 paired rows:  {len(v7)}  (Functional='No')")
    print(f"v49 paired rows: {len(v49)} (Functional='Yes')")

    combined = pd.concat([v7, v49], ignore_index=True)
    print(f"Combined: {len(combined)} rows x {len(combined.columns)} cols")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    combined.to_csv(args.output, index=False)
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
