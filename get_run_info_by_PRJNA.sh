#!/bin/bash
#
# Get run-info and metadata for a single PRJNA ID
#
# -----------------------------------------------------------------------------
set -Eeuo pipefail

BRIOPROJECT=$1
OUT_DIR=$2

# get the SRA run-info for a single PROJECT ID
esearch -db sra -query \"$BIOPROJECT\" | efetch -format runinfo > $OUT_DIR/run-info.csv

# Get additional metadata not present in run_info.csv
esearch -db sra -query \"$BIOPROJECT\" | \
    efetch -format native | \
    xtract -pattern EXPERIMENT_PACKAGE \
        -group RUN_SET \
            -block RUN \
                -block IDENTIFIERS \
                    -element PRIMARY_ID \
        -group SAMPLE \
            -block SAMPLE_ATTRIBUTES \
                -block SAMPLE_ATTRIBUTE \
                    -element VALUE \
        -group EXPERIMENT \
            -element TITLE > $OUT_DIR/metadata.txt
