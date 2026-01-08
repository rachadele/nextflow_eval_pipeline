#!/bin/bash
set -e

# Run human test pipeline
nextflow run main.nf \
    -profile conda \
    -params-file params.hs.json \
    -work-dir hsap \
    -resume
