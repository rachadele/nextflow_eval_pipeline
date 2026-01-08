#!/bin/bash
set -e

# Run mouse test pipeline
nextflow run main.nf \
    -profile conda \
    -params-file params.mm.json \
    -work-dir mmus \
    -resume 
