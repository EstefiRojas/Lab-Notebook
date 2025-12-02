#!/bin/bash

# Stop on errors
set -e

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_csv_file> <gencode_gtf_file>"
    exit 1
fi

INPUT_FILE=$1
GTF_FILE=$2
TEMP_PRELIM="temp_preliminary.tsv"
TEMP_MISSING="temp_missing_ids.txt"
TEMP_API_MAP="temp_api_mapping.tsv"
FINAL_OUTPUT="../results/annotated_unified_genome_alignments.csv"

echo "========================================================"
echo "   Gene Annotator: GTF + Coordinate-Based API Lookup"
echo "========================================================"

# Determine unzipping method
if [[ "$GTF_FILE" == *.gz ]]; then CAT_CMD="gunzip -c"; else CAT_CMD="cat"; fi

# ---------------------------------------------------------
# PHASE 1: Local GTF Mapping & Coordinate Extraction
# ---------------------------------------------------------
echo "[1/4] Mapping names and extracting coordinates from GTF..."

# We use AWK to build two maps:
# 1. ID -> Name
# 2. ID -> Chr:Start:End
awk '
BEGIN { OFS="," }
FNR==NR {
    if ($3 == "gene") {
        gid=""; gname=""
        
        # Extract attributes
        if (match($9, /gene_id "[^"]+"/)) gid = substr($9, RSTART+9, RLENGTH-10)
        if (match($9, /gene_name "[^"]+"/)) gname = substr($9, RSTART+11, RLENGTH-12)
        
        if (gid != "") {
            split(gid, a, "."); base=a[1]
            
            # Save Name
            if (gname != "") map[base] = gname
            
            # Save Coordinates (Chr, Start, End)
            # $1=Chr, $4=Start, $5=End
            coords[base] = $1 "\t" $4 "\t" $5
        }
    }
    next
}
{
    if (FNR == 1) { 
        sub(/\r$/, "", $0)
        print $0, "Gene_Name" 
    }
    else {
        sub(/\r$/, "", $0)
        # CSV Input: ENSG_ID is in column 6
        split($6, a, "."); id=a[1]
        val = "NA"
        
        # Try to name it locally
        if (id in map) {
            if (map[id] != id && map[id] !~ /^ENSG/) {
                val = map[id]
            }
        }
        print $0, val
    }
    
    # Side channel: Create a coordinate file for missing items
    # We only write to this if we failed to find a name
    if (FNR > 1) {
        split($6, a, "."); id=a[1]
        val = "NA"
        if (id in map && map[id] != id && map[id] !~ /^ENSG/) val = map[id]
        
        if (val == "NA" && (id in coords)) {
            # Write to specific file: ID, Chr, Start, End
            # Explicitly use tabs for Python script compatibility
            print id "\t" coords[id] > "'"$TEMP_MISSING"'"
        }
    }
}
' FS='\t' <($CAT_CMD "$GTF_FILE") FS=',' "$INPUT_FILE" > "$TEMP_PRELIM"

# ---------------------------------------------------------
# PHASE 2: Check Missing
# ---------------------------------------------------------
echo "[2/4] identifying missing gene names..."

if [ ! -f "$TEMP_MISSING" ] || [ ! -s "$TEMP_MISSING" ]; then
    echo "   > Success: All genes named locally."
    mv "$TEMP_PRELIM" "$FINAL_OUTPUT"
    exit 0
fi

# Deduplicate the missing file
sort -u "$TEMP_MISSING" -o "$TEMP_MISSING"
MISSING_COUNT=$(wc -l < "$TEMP_MISSING")
echo "   > $MISSING_COUNT IDs require API lookup."

# ---------------------------------------------------------
# PHASE 3: Python Multi-API Waterfall (Coords Supported)
# ---------------------------------------------------------
echo "[3/4] Querying APIs..."

python3 - "$TEMP_MISSING" "$TEMP_API_MAP" << 'EOF'
import sys
import json
import time
import urllib.request
import urllib.parse
import ssl

# Bypass SSL
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

infile = sys.argv[1]
outfile = sys.argv[2]

# Load Missing IDs and Coordinates into a dictionary
# Format: { "ENSG...": { "chr": "1", "start": "100", "end": "200" } }
missing_data = {}
with open(infile) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 4:
            eid = parts[0]
            missing_data[eid] = {
                "chr": parts[1],
                "start": parts[2],
                "end": parts[3]
            }

results = {} 

def print_progress(current, total, prefix):
    percent = (current / total) * 100
    sys.stderr.write(f"\r     {prefix}: {current}/{total} ({percent:.1f}%)")
    sys.stderr.flush()

def get_remaining_ids():
    return [x for x in missing_data.keys() if x not in results]

# ==============================================================================
# 1. MyGene.info (Batch)
# ==============================================================================
remaining = get_remaining_ids()
if remaining:
    sys.stderr.write(f"\n   > MyGene.info (Batch): {len(remaining)} IDs\n")
    mg_url = "https://mygene.info/v3/query"
    batch_size = 1000
    
    # Convert keys to list for batching
    rem_list = list(remaining)
    total = len(rem_list)
    
    for i in range(0, total, batch_size):
        print_progress(i, total, "Batching")
        batch = rem_list[i:i+batch_size]
        data = {
            "q": ",".join(batch),
            "scopes": "ensembl.gene,alias",
            "fields": "symbol,name,alias",
            "species": "human"
        }
        try:
            encoded_data = urllib.parse.urlencode(data).encode("utf-8")
            req = urllib.request.Request(mg_url, data=encoded_data)
            with urllib.request.urlopen(req, context=ctx) as res:
                hits = json.loads(res.read().decode())
                for hit in hits:
                    if "symbol" in hit and not hit["symbol"].startswith("ENSG"):
                        results[hit["query"]] = hit["symbol"]
        except Exception:
            pass
    print_progress(total, total, "Batching")

# ==============================================================================
# 2. Ensembl REST API (ID Lookup)
# ==============================================================================
#remaining = get_remaining_ids()
#if remaining:
#    sys.stderr.write(f"\n   > Ensembl API: {len(remaining)} IDs\n")
#    total = len(remaining)
#    
#    for i, ens_id in enumerate(remaining, 1):
#        print_progress(i, total, "Querying")
#        try:
#            url = f"https://rest.ensembl.org/lookup/id/{ens_id}?content-type=application/json"
#            req = urllib.request.Request(url)
#            with urllib.request.urlopen(req, context=ctx) as res:
#                data = json.loads(res.read().decode())
#                if "display_name" in data and data["display_name"]:
#                     name = data["display_name"]
#                     if not name.startswith("ENSG"):
#                        results[ens_id] = name
#            time.sleep(0.07) 
#        except Exception:
#            pass
#
# ==============================================================================
# 3. LncBook (Search by Genomic Position)
# ==============================================================================
#remaining = get_remaining_ids()
#if remaining:
#    sys.stderr.write(f"\n   > LncBook (Positional Search): {len(remaining)} IDs\n")
#    total = len(remaining)
#    
#    for i, ens_id in enumerate(remaining, 1):
#        print_progress(i, total, "Searching")
#        
#        coords = missing_data[ens_id]
#        chrom = coords['chr']
#        
#        # LncBook requires 'chr' prefix (e.g., chr1, chrX)
#        if not chrom.startswith("chr"):
#            chrom = f"chr{chrom}"
#            
#        start = coords['start']
#        end = coords['end']
#
#        try:
#            # Endpoint: Query genes by Region
#            # Note: We use the LncBook API specific to Genes, not EWAS
#            url = f"https://ngdc.cncb.ac.cn/lncbook/api/gene/region?chromosome={chrom}&start={start}&end={end}"
#           
#            req = urllib.request.Request(url)
#            with urllib.request.urlopen(req, context=ctx) as res:
#                data = json.loads(res.read().decode())
#                
#                # Structure: { "data": [ { "symbol": "LINC00...", ... } ], ... }
#                # Or direct list depending on version. Based on current NGDC API:
#                if "data" in data and len(data["data"]) > 0:
#                    # We take the first result that has a symbol
#                    for gene in data["data"]:
#                        if "symbol" in gene and gene["symbol"]:
#                            # Optional: Check overlap strictly? 
#                           # For now, assume the region query implies overlap
#                            results[ens_id] = gene["symbol"]
#                            break
#                            
#            time.sleep(0.2)
#        except Exception:
#            pass
#
sys.stderr.write("\n")

# Output
with open(outfile, "w") as out:
    for eid, name in results.items():
        out.write(f"{eid}\t{name}\n")
EOF

# ---------------------------------------------------------
# PHASE 4: Merging
# ---------------------------------------------------------
echo -e "\n[4/4] Merging annotations..."

if [ -s "$TEMP_API_MAP" ]; then
    awk '
    BEGIN { OFS="," }
    FNR==NR { map[$1] = $2; next }
    {
        if (FNR == 1) { 
            sub(/\r$/, "", $0)
            print $0; next 
        }
        sub(/\r$/, "", $0)
        if ($NF == "NA") {
            split($6, a, "."); id=a[1]
            if (id in map) $NF = map[id]
        }
        print $0
    }
    ' FS='\t' "$TEMP_API_MAP" FS=',' "$TEMP_PRELIM" > "$FINAL_OUTPUT".tmp
    
    # Reorder columns
    # Current (Input to this step): 
    # 1:Study, 2:gRNA_ID, 3:Target_Gene_ID, 4:gRNA_Type, 5:Essentiality, 
    # 6:ENSG_ID, 7:Strand, 8:Probability_Functional, 9:Protein_Off_Target, 
    # 10:gRNA_Sequence, 11:Antisense_to_CDS, 12:Antisense_Target_Name, 13:Gene_Name
    
    # Desired Order:
    # 1:Study
    # 2:gRNA_ID
    # 10:gRNA_Sequence
    # 4:gRNA_Type
    # 3:Target_Gene_ID (Original Target)
    # 6:ENSG_ID (Mapped Target)
    # 13:Gene_Name (Mapped Name)
    # 7:Strand
    # 5:Essentiality
    # 8:Probability_Functional
    # 9:Protein_Off_Target
    # 11:Antisense_to_CDS
    # 12:Antisense_Target_Name
    
    awk -F',' 'BEGIN{OFS=","} {
        # Handle header
        if (NR==1) {
            print "Study", "gRNA_ID", "gRNA_Sequence", "gRNA_Type", "Target_Gene_ID", "ENSG_ID", "Gene_Name", "Strand", "Essentiality", "Probability_Functional", "Protein_Off_Target", "Antisense_to_CDS", "Antisense_Target_Name"
        } else {
            print $1, $2, $10, $4, $3, $6, $13, $7, $5, $8, $9, $11, $12
        }
    }' "$FINAL_OUTPUT".tmp > "$FINAL_OUTPUT"
    rm "$FINAL_OUTPUT".tmp
else
    echo "   > No new names found. Keeping NAs."
    awk -F',' 'BEGIN{OFS=","} { 
        sub(/\r$/, "", $0)
        # Reorder even if no new names (Gene_Name is still added as NA/Gene_Name in previous step)
        if (NR==1) {
             print "Study", "gRNA_ID", "gRNA_Sequence", "gRNA_Type", "Target_Gene_ID", "ENSG_ID", "Gene_Name", "Strand", "Essentiality", "Probability_Functional", "Protein_Off_Target", "Antisense_to_CDS", "Antisense_Target_Name"
        } else {
             print $1, $2, $10, $4, $3, $6, $13, $7, $5, $8, $9, $11, $12
        }
    }' "$TEMP_PRELIM" > "$FINAL_OUTPUT"
fi

rm -f "$TEMP_PRELIM" "$TEMP_MISSING" "$TEMP_API_MAP"
echo "Done. Annotated file saved to: $FINAL_OUTPUT"