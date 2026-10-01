#!/bin/bash
# Test RF + kNN classifiers on both organisms with a small fixed subsample.
set -e

SUBSAMPLE_REF=100
SUBSAMPLE_QUERY=100
KNN_K=15
OUTDIR_PREFIX="test_knn"

echo "=== Running human (homo_sapiens) ==="
nextflow run main.nf \
    -profile conda \
    -params-file params.hs.json \
    -work-dir work_test_knn_hsap \
    --subsample_ref ${SUBSAMPLE_REF} \
    --subsample_query ${SUBSAMPLE_QUERY} \
    --knn_n_neighbors ${KNN_K} \
    --outdir_prefix "${OUTDIR_PREFIX}/hsap" \
    -resume

echo "=== Running mouse (mus_musculus) ==="
nextflow run main.nf \
    -profile conda \
    -params-file params.mm.json \
    -work-dir work_test_knn_mmus \
    --subsample_ref ${SUBSAMPLE_REF} \
    --subsample_query ${SUBSAMPLE_QUERY} \
    --knn_n_neighbors ${KNN_K} \
    --outdir_prefix "${OUTDIR_PREFIX}/mmus" \
    -resume

echo "=== Done. Results in ${OUTDIR_PREFIX}/ ==="
