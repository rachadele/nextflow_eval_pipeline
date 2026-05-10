#!/bin/bash
# Run on branch: exclude-ref-unsupported-metrics
set -e

nextflow run main.nf \
    -profile conda \
    -params-file params.hs.json \
    -work-dir work_excl \
    --outdir_prefix test_excl \
    -resume
