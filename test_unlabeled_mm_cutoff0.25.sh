#!/bin/bash
# Test: mouse, cutoff 0.25, queries augmented with author-unlabeled cells.
# Builds the augmented h5ads (skipped if already built), then runs one nextflow job.
# Everything goes to cosmos. Run from the worktree root, inside tmux.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

# The build keeps all unlabeled cells per sample (native ratios). Ratios per study: /space/grp/rschwartz/rschwartz/mask-annotation-overlap/results/qc_labeled_vs_unlabeled/unlabeled_ratio_by_study.tsv (scripts/qc_labeled_vs_unlabeled.py, commit c859220)
cutoff=0.25
subsample_ref=500
subsample_query=100
ref_split=dataset_id
normalization_method=SCT
use_gap=false

base=/cosmos/data/nextflow-eval-pipeline
queries_dir=$base/unlabeled_queries/mouse_ratio_native
outdir=$base/results/census-map-fixes_with_unlabeled_test/mus_musculus/ratio_native/ref_$subsample_ref/cutoff_$cutoff

if [ ! -d "$queries_dir" ]; then
    /home/rschwartz/anaconda3/envs/scanpyenv/bin/python scripts/build_unlabeled_queries.py --dest "$queries_dir"
fi

nextflow main.nf -params-file params.mm.json \
    --queries_adata "$queries_dir/h5ad/**h5ad" \
    --outdir "$outdir" \
    --subsample_query "$subsample_query" \
    --subsample_ref "$subsample_ref" \
    --ref_split "$ref_split" \
    -profile conda \
    --cutoff "$cutoff" \
    --subset_type sample \
    --batch_correct true \
    -resume \
    --remove_unknown true \
    --normalization_method "$normalization_method" \
    -process.executor slurm \
    --use_gap "$use_gap"

find "$outdir" -name '*.unlabeled_unknown_rate.*.tsv.gz'
