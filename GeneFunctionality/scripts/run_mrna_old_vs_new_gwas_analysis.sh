#!/usr/bin/env bash
#
# End-to-end GWAS SNP density analysis for the mrna_old_vs_new paired dataset
# (250 v7-deprecated genes + 250 v49-current genes, each with exon2 + exon3).
#
# Steps:
#   1. Build wide-format CSV from predicting-genic-features predictions.
#   2. Resolve transcript coordinates per source (v7 GTF for v7 rows,
#      v49 GTF for v49 rows), then concatenate.
#   3. Intersect with GWAS catalog via extract_gwas_data_mrna.sh.
#   4. Run R analysis to produce plots.
#
# Run from the GeneFunctionality directory:
#   bash scripts/run_mrna_old_vs_new_gwas_analysis.sh

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate base

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_DIR"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PGF="${PROJECT_DIR}/../predicting-genic-features"

WIDE_NO_TR="data/model_predictions/mrna_old_vs_new_exon2-exon3-features.csv"
WIDE_WITH_TR="data/model_predictions/mrna_old_vs_new_exon2-exon3-features-with-transcripts.csv"
GWAS_OUT="results/mrna_old_vs_new_gwas_snp_density.csv"

GTF_V49="data/references/gencode.v49.primary_assembly.annotation.gtf"
# v7 GTF is GRCh37; the dataset uses hg38 coordinates. Use the pre-lifted hg38
# protein-coding gene BED instead, so all 250 v7 rows resolve. (Falls back to
# raw GRCh37 GTF if the lifted BED is unavailable, which will leave ~113/250
# v7 rows with NA transcript coords due to assembly mismatch.)
V7_HG38_BED="${PGF}/data/datasets/mrna_old_vs_new/_tmp_v7_deprecated/v7_pc_genes.hg38.sorted.bed"
GTF_V7="${PGF}/data/datasets/gencode_v7_mrna/gencode.v7.annotation.gtf"
GWAS_CATALOG="data/gwas/gwas_catalog_v1.0-associations_e114_r2025-07-21.tsv"

# ---------------------------------------------------------------------------
# Step 1: Build wide CSV
# ---------------------------------------------------------------------------
echo "[Step 1/4] Building wide CSV from predicting-genic-features predictions..."
python3 scripts/build_mrna_old_vs_new_wide.py \
    --features-dir "${PGF}/results/mrna_old_vs_new" \
    --output "${WIDE_NO_TR}"

# ---------------------------------------------------------------------------
# Step 2: Resolve transcript coordinates per source (v7 -> v7 GTF, v49 -> v49 GTF)
# ---------------------------------------------------------------------------
echo ""
echo "[Step 2/4] Resolving transcript coordinates per source..."

WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

V7_IN="${WORK_DIR}/v7_only.csv"
V49_IN="${WORK_DIR}/v49_only.csv"
V7_OUT="${WORK_DIR}/v7_with_transcripts.csv"
V49_OUT="${WORK_DIR}/v49_with_transcripts.csv"

# Split by ID prefix
head -1 "$WIDE_NO_TR" > "$V7_IN"
awk -F',' 'NR==1 {next} $1 ~ /^v7_gene_/' "$WIDE_NO_TR" >> "$V7_IN"

head -1 "$WIDE_NO_TR" > "$V49_IN"
awk -F',' 'NR==1 {next} $1 ~ /^v49_gene_/' "$WIDE_NO_TR" >> "$V49_IN"

echo "  v7 rows:  $(( $(wc -l < "$V7_IN") - 1 ))"
echo "  v49 rows: $(( $(wc -l < "$V49_IN") - 1 ))"

if [ -f "$V7_HG38_BED" ]; then
    echo "  Running get_mrna_transcript_coords.py on v7 against hg38-lifted v7 BED..."
    python3 scripts/get_mrna_transcript_coords.py \
        --input "$V7_IN" \
        --bed "$V7_HG38_BED" \
        --output "$V7_OUT"
else
    echo "  WARNING: hg38-lifted v7 BED missing — falling back to GRCh37 GTF (expect ~113 unresolved)"
    python3 scripts/get_mrna_transcript_coords.py \
        --input "$V7_IN" \
        --gtf "$GTF_V7" \
        --output "$V7_OUT"
fi

echo "  Running get_mrna_transcript_coords.py on v49 against gencode.v49..."
python3 scripts/get_mrna_transcript_coords.py \
    --input "$V49_IN" \
    --gtf "$GTF_V49" \
    --output "$V49_OUT"

# Concatenate (header from v7 file, append v49 data rows)
mkdir -p "$(dirname "$WIDE_WITH_TR")"
cp "$V7_OUT" "$WIDE_WITH_TR"
tail -n +2 "$V49_OUT" >> "$WIDE_WITH_TR"

echo "  Combined: $(( $(wc -l < "$WIDE_WITH_TR") - 1 )) rows -> $WIDE_WITH_TR"

# Quick sanity: how many v7 rows got a transcript span (i.e. resolved gene)?
V7_RESOLVED=$(awk -F',' 'NR>1 && $1 ~ /^v7_gene_/ && $40 != "NA" {c++} END {print c+0}' "$WIDE_WITH_TR")
V49_RESOLVED=$(awk -F',' 'NR>1 && $1 ~ /^v49_gene_/ && $40 != "NA" {c++} END {print c+0}' "$WIDE_WITH_TR")
echo "  Transcript-coord resolution: v7=${V7_RESOLVED}/250  v49=${V49_RESOLVED}/250"

# ---------------------------------------------------------------------------
# Step 3: GWAS intersect
# ---------------------------------------------------------------------------
echo ""
echo "[Step 3/4] Running extract_gwas_data_mrna.sh..."
bash scripts/extract_gwas_data_mrna.sh "$WIDE_WITH_TR" "$GWAS_CATALOG" "$GWAS_OUT"

# ---------------------------------------------------------------------------
# Step 3.5: Append the 100 v49 positive-control mRNAs (already in GWAS-output form)
# into the v49 group of the final density CSV so all downstream plots include them.
# ---------------------------------------------------------------------------
EXTRA_100="results/mrna_100_gwas_snp_density.csv"
if [ -f "$EXTRA_100" ]; then
    echo ""
    echo "[Step 3.5/4] Appending $EXTRA_100 (100 v49 controls) into final density CSV..."
    python3 - <<PYEOF
import pandas as pd

main_path  = "$GWAS_OUT"
extra_path = "$EXTRA_100"

main  = pd.read_csv(main_path)
extra = pd.read_csv(extra_path)
print(f"  Main:  {len(main)} rows x {len(main.columns)} cols")
print(f"  Extra: {len(extra)} rows x {len(extra.columns)} cols")

# Align by column name (the two files have the same names but different orders);
# pandas concat will reorder extra to match main automatically.
combined = pd.concat([main, extra[main.columns.tolist()]], ignore_index=True)
print(f"  Combined: {len(combined)} rows  (Functional split: {combined['Functional'].value_counts().to_dict()})")
combined.to_csv(main_path, index=False)
print(f"  Wrote: {main_path}")
PYEOF
else
    echo "  WARNING: $EXTRA_100 not found — skipping 100-control append"
fi

# ---------------------------------------------------------------------------
# Step 4: R analysis
# ---------------------------------------------------------------------------
echo ""
echo "[Step 4/4] Running R analysis..."
mkdir -p results/gwas
cd scripts
Rscript mrna_old_vs_new_gwas_snp_density_analysis.R
cd "$PROJECT_DIR"

echo ""
echo "=========================================="
echo " Done."
echo "=========================================="
echo "GWAS density CSV: $GWAS_OUT"
echo "Plots:"
ls -1 results/gwas/mrna_old_vs_new_*.png 2>/dev/null
