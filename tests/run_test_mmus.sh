#!/bin/bash
set -e

# Run mouse test pipeline
nextflow run main.nf \
    -profile conda,test_mmus \
    -work-dir mmus \
    -resume \
    "$@"
