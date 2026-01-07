#!/bin/bash
set -e

# Run human test pipeline
nextflow run main.nf \
    -profile conda,test_hsap \
    -work-dir hsap \
    -resume \
    "$@"
