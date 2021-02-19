#!/bin/bash
#
# Query the SRA database for each Project Accession returned by the search
# script. Return the run-info.csv file for each which contains information
# about the SRA Accessions, treatments, library prep, etc.
#
# -----------------------------------------------------------------------------
set -Eeuo pipefail

SEARCH_RESULTS=$1
OUT_DIR=$2
MIN_SAMPLES=1
MAX_SAMPLES=100

while read -r BIOPROJECT N_SAMPLES TITLE; do
    if [ $N_SAMPLES -gt $MIN_SAMPLES -a $N_SAMPLES -lt $MAX_SAMPLES ] 
    then
        mkdir -p $OUT_DIR/$BIOPROJECT
        sleep 1
        esearch -db sra -query \"$BIOPROJECT\" < /dev/null | efetch -format runinfo > $OUT_DIR/$BIOPROJECT/run-info.csv
    fi
done < $SEARCH_RESULTS
