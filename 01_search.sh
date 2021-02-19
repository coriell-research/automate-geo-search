#!/bin/bash
#
# Search NCBI for cancer datasets. Return a flat file with the BioProject
# Accession, Number of Samples, and Title of Project
# 
# QUERY below can eventually be replace with a cmdline arg
# -----------------------------------------------------------------------------
set -Eeuo pipefail

QUERY=$1
OUT_DIR=$2

esearch \
    -db gds \
    -query "$QUERY" | \
efetch \
    -format docsum | \
xtract \
    -pattern DocumentSummary \
    -element BioProject n_samples title > $OUT_DIR/search-results.tsv