#!/bin/bash
# Run on branch: main
set -e

nextflow run main.nf \
    -profile conda \
    -params-file params.hs.json \
    -work-dir work_main \
    --outdir_prefix test_main \
    -resume
