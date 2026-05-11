#!/bin/bash
# Run on branch: ma-et-al-support
set -e

nextflow run main.nf \
    -profile conda \
    -params-file params.hs.ma.json \
    -work-dir work_ma \
    -resume
