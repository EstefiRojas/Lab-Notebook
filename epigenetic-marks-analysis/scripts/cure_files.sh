#!/bin/bash
set -uex

CURE_FILES_FOLDER=$1

awk -F'\t' '$5 == "fold change over control" && $6 == "GRCh38" && $56 == "released" {print $48}' "$CURE_FILES_FOLDER"/metadata.tsv > "$CURE_FILES_FOLDER"/files.txt
