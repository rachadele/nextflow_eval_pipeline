#!/bin/bash

nextflow run main.nf \
    -profile conda \
    -params-file params.hs.json \
    --subsample_query 100 \
    --outdir_prefix "results/hsap/test_ma_et_al" \
    -resume
