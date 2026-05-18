#!/usr/bin/env python3
"""
Get transcript-level coordinates for mRNA genes from GENCODE v49 GTF.

Uses bedtools intersect to overlap exon coordinates with GENCODE gene features,
then extracts the parent gene boundaries as transcript coordinates.

Usage:
    python get_mrna_transcript_coords.py \
        --input data/model_predictions/matrix-positive-100mrna_exon2-exon3-features.csv \
        --gtf data/references/gencode.v49.primary_assembly.annotation.gtf \
        --output data/model_predictions/matrix-positive-100mrna_exon2-exon3-features-with-transcripts.csv
"""

import argparse
import csv
import subprocess
import tempfile
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add transcript coordinates from GENCODE GTF to mRNA feature file."
    )
    parser.add_argument("--input", required=True, help="Input mRNA CSV file")
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--gtf", help="GENCODE GTF annotation file")
    src.add_argument("--bed",
                     help="Pre-built sorted gene BED (chr, start, end, gene_name[|gene_id], score, strand). "
                          "Used as-is; assumes coordinates are already on the target assembly.")
    parser.add_argument("--output", required=True, help="Output CSV file with transcript columns")
    return parser.parse_args()


def adopt_bed_genes(bed_path, tmpdir):
    """Normalize a pre-built gene BED to the format extract_gtf_genes produces.

    Accepts column 4 in either form: 'gene_name' or 'gene_name|gene_id' (the
    latter is the convention used by v7_pc_genes.hg38.sorted.bed). Re-sorts to
    be safe.
    """
    out_path = os.path.join(tmpdir, "gencode_genes.bed")
    count = 0
    with open(bed_path, "r") as bed, open(out_path, "w") as out:
        for line in bed:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end = parts[0], parts[1], parts[2]
            name_field = parts[3]
            # Strip a '|ENSG...' suffix if present
            gene_name = name_field.split("|")[0]
            strand = parts[5] if len(parts) > 5 else "."
            out.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, start, end, gene_name, strand))
            count += 1

    sorted_path = os.path.join(tmpdir, "gencode_genes_sorted.bed")
    subprocess.run(
        "sort -k1,1 -k2,2n {} > {}".format(out_path, sorted_path),
        shell=True, check=True
    )
    print("  -> Adopted {} gene entries from pre-built BED".format(count))
    return sorted_path


def extract_gtf_genes(gtf_path, tmpdir):
    """Parse GENCODE GTF for protein_coding gene features, write as sorted BED."""
    bed_path = os.path.join(tmpdir, "gencode_genes.bed")
    count = 0

    with open(gtf_path, "r") as gtf, open(bed_path, "w") as bed:
        for line in gtf:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "gene":
                continue

            attributes = parts[8]
            if 'gene_type "protein_coding"' not in attributes:
                continue

            chrom = parts[0]
            start = int(parts[3]) - 1  # Convert to 0-based
            end = int(parts[4])
            strand = parts[6]

            # Extract gene_name
            gene_name = "Unknown"
            if 'gene_name "' in attributes:
                gene_name = attributes.split('gene_name "')[1].split('"')[0]

            bed.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, start, end, gene_name, strand))
            count += 1

    # Sort the BED file
    sorted_path = os.path.join(tmpdir, "gencode_genes_sorted.bed")
    subprocess.run(
        "sort -k1,1 -k2,2n {} > {}".format(bed_path, sorted_path),
        shell=True, check=True
    )
    print("  -> Extracted {} protein-coding genes from GENCODE GTF".format(count))
    return sorted_path


def create_exon_bed(input_csv, tmpdir):
    """Create BED file from exon2 coordinates for intersection."""
    bed_path = os.path.join(tmpdir, "exons.bed")
    genes = []

    with open(input_csv, "r") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            chrom = row["ex2_Chromosome"].strip()
            start = int(row["ex2_Start"].strip())
            end = int(row["ex2_End"].strip())

            # Normalize chromosome
            if not chrom.startswith("chr"):
                chrom = "chr{}".format(chrom)

            # Store gene info for later mapping
            genes.append({
                "line_num": i + 1,
                "chrom": chrom,
                "start": start,
                "end": end,
                "gene_id": row.get("GeneID", "").strip(),
            })

    # Write BED (0-based start)
    with open(bed_path, "w") as bed:
        for g in genes:
            bed.write("{}\t{}\t{}\t{}\n".format(g['chrom'], g['start'] - 1, g['end'], g['line_num']))

    sorted_path = os.path.join(tmpdir, "exons_sorted.bed")
    subprocess.run(
        "sort -k1,1 -k2,2n {} > {}".format(bed_path, sorted_path),
        shell=True, check=True
    )
    print("  -> Created BED file for {} exon2 regions".format(len(genes)))
    return sorted_path, genes


def run_intersect(exon_bed, gene_bed, tmpdir):
    """Run bedtools intersect to find overlapping genes."""
    output_path = os.path.join(tmpdir, "intersect.bed")
    cmd = "bedtools intersect -wa -wb -a {} -b {} > {}".format(exon_bed, gene_bed, output_path)
    result = subprocess.run(cmd, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode != 0:
        print("  -> bedtools error: {}".format(result.stderr.decode()), file=sys.stderr)
        sys.exit(1)

    # Parse results: map line_num -> gene info
    transcript_map = {}
    with open(output_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            line_num = int(parts[3])
            gene_chrom = parts[4]
            gene_start = int(parts[5]) + 1  # Convert back to 1-based
            gene_end = int(parts[6])
            gene_name = parts[7]

            # If multiple genes overlap, keep the one with the largest overlap
            if line_num not in transcript_map:
                transcript_map[line_num] = {
                    "chrom": gene_chrom,
                    "start": gene_start,
                    "end": gene_end,
                    "gene_name": gene_name,
                }
            else:
                # Keep the gene with the larger span (more likely to be the actual parent)
                existing = transcript_map[line_num]
                existing_span = existing["end"] - existing["start"]
                new_span = gene_end - gene_start
                if new_span > existing_span:
                    transcript_map[line_num] = {
                        "chrom": gene_chrom,
                        "start": gene_start,
                        "end": gene_end,
                        "gene_name": gene_name,
                    }

    matched = len(transcript_map)
    print("  -> Matched {} genes via bedtools intersect".format(matched))
    return transcript_map


def write_output(input_csv, output_csv, transcript_map):
    """Write output CSV with transcript coordinate columns appended."""
    with open(input_csv, "r") as fin, open(output_csv, "w", newline="") as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout)

        # Write header
        header = next(reader)
        header.extend([
            "Chrm_Transcript", "Start_Transcript", "End_Transcript", "gene_name"
        ])
        writer.writerow(header)

        # Write data rows
        missing = 0
        for i, row in enumerate(reader):
            line_num = i + 1
            if line_num in transcript_map:
                t = transcript_map[line_num]
                row.extend([t["chrom"], str(t["start"]), str(t["end"]), t["gene_name"]])
            else:
                row.extend(["NA", "NA", "NA", "NA"])
                missing += 1

            writer.writerow(row)

        if missing > 0:
            print("  -> WARNING: {} genes had no GENCODE match".format(missing))


def main():
    args = parse_args()

    if not os.path.exists(args.input):
        print("Error: Input file '{}' not found.".format(args.input), file=sys.stderr)
        sys.exit(1)
    gene_src = args.gtf if args.gtf else args.bed
    if not os.path.exists(gene_src):
        print("Error: gene-source file '{}' not found.".format(gene_src), file=sys.stderr)
        sys.exit(1)

    # Check bedtools
    result = subprocess.run("command -v bedtools", shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print("Error: bedtools is required but not installed.", file=sys.stderr)
        sys.exit(1)

    print("=" * 50)
    print("Get mRNA Transcript Coordinates from GENCODE")
    print("=" * 50)

    with tempfile.TemporaryDirectory() as tmpdir:
        if args.gtf:
            print("\nStep 1: Parsing GENCODE GTF for protein-coding genes...")
            gene_bed = extract_gtf_genes(args.gtf, tmpdir)
        else:
            print("\nStep 1: Adopting pre-built gene BED ({})...".format(args.bed))
            gene_bed = adopt_bed_genes(args.bed, tmpdir)

        print("\nStep 2: Creating BED file from exon2 coordinates...")
        exon_bed, genes = create_exon_bed(args.input, tmpdir)

        print("\nStep 3: Running bedtools intersect...")
        transcript_map = run_intersect(exon_bed, gene_bed, tmpdir)

        print("\nStep 4: Writing output with transcript columns...")
        write_output(args.input, args.output, transcript_map)

    print("\nDone. Output written to: {}".format(args.output))


if __name__ == "__main__":
    main()
