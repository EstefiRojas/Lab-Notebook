#!/usr/bin/env bash
# Add (+) or (-) to the end of each FASTA header using strand from a GENCODE GTF.
# Works with BSD awk / POSIX awk (no gawk-only features).
# Usage:
#   ./add_strand_to_fasta.sh [-v] [-a] input.fa[.gz] annotations.gtf[.gz] > output.fa
set -euo pipefail

VERBOSE=0
ALWAYS=0
while getopts ":va" opt; do
  case "$opt" in
    v) VERBOSE=1 ;;
    a) ALWAYS=1 ;;
    *) echo "Usage: $0 [-v] [-a] <input.fasta[.gz]> <annotations.gtf[.gz]>" >&2; exit 1 ;;
  esac
done
shift $((OPTIND-1))

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 [-v] [-a] <input.fasta[.gz]> <annotations.gtf[.gz]>" >&2
  exit 1
fi

FASTA="$1"
GTF="$2"

read_stream() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    gzip -cd -- "$f"
  else
    cat -- "$f"
  fi
}

awk -v VERBOSE="$VERBOSE" -v ALWAYS="$ALWAYS" '
# ---------- small helpers ----------
function trim(s){ gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", s); return s }

# Parse GTF attributes into global array attrs[key]=value (POSIX-safe)
# Accepts chunks like: key "value"
function parse_gtf_attrs(attr_str,   n,i,chunk,sp,firstq,rest,secondq,key,val,parts) {
  delete attrs
  n = split(attr_str, parts, ";")
  for (i=1; i<=n; i++) {
    chunk = trim(parts[i])
    if (chunk == "") continue
    # must contain a space then a quoted value
    sp = index(chunk, " ")
    firstq = index(chunk, "\"")
    if (sp == 0 || firstq == 0) continue
    key = substr(chunk, 1, sp-1)
    rest = substr(chunk, firstq+1)           # after first quote
    secondq = index(rest, "\"")
    if (secondq == 0) continue
    val = substr(rest, 1, secondq-1)
    attrs[key] = val
  }
}

# ---------- PASS 1: GTF -> map transcript_id(.version) => strand ----------
FNR==NR {
  if ($0 ~ /^#/) next
  nf = split($0, col, "\t")
  if (nf < 9) next

  feature = col[3]
  strand  = col[7]
  if (strand != "+" && strand != "-") next

  attrs_raw = col[9]
  parse_gtf_attrs(attrs_raw)

  tid = ("transcript_id" in attrs) ? attrs["transcript_id"] : ""
  if (tid == "") next
  tver = ("transcript_version" in attrs) ? attrs["transcript_version"] : ""

  if (feature == "transcript" || !(tid in t2s)) {
    t2s[tid] = strand
    if (tver != "") t2s[tid "." tver] = strand
  }
  next
}

# ---------- PASS 2: FASTA -> append (+) or (-) ----------
{
  if (substr($0,1,1) == ">") {
    header = substr($0, 2)
    firstField = header
    p = index(header, "|")
    if (p > 0) firstField = substr(header, 1, p-1)

    tid = firstField
    strand = ""
    if (tid in t2s) {
      strand = t2s[tid]
    } else {
      dotp = index(tid, ".")
      if (dotp > 0) {
        tid_nover = substr(tid, 1, dotp-1)
        if (tid_nover in t2s) strand = t2s[tid_nover]
      }
    }

    if (strand != "") {
      sep = (header ~ /\|$/) ? "" : "|"
      print ">" header sep "(" strand ")"
      changed++
    } else if (ALWAYS) {
      sep = (header ~ /\|$/) ? "" : "|"
      print ">" header sep "(?)"
      changed_unknown++
    } else {
      print $0
      unchanged++
    }
    total_headers++
  } else {
    print $0
  }
}
END{
  if (VERBOSE) {
    k=0; for (x in t2s) k++
    print "Loaded transcript->strand keys: " k > "/dev/stderr"
    print "FASTA headers seen: " (total_headers+0) > "/dev/stderr"
    print "Headers updated (+/-): " (changed+0) > "/dev/stderr"
    if (ALWAYS) print "Headers marked (?): " (changed_unknown+0) > "/dev/stderr"
    print "Headers left unchanged: " (unchanged+0) > "/dev/stderr"
  }
}
' <(read_stream "$GTF") <(read_stream "$FASTA")
