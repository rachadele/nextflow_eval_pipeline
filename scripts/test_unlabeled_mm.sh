#!/bin/bash
# Test: mouse, queries augmented with author-unlabeled cells. Usage: scripts/test_unlabeled_mm.sh [cutoff] (default 0.25).
# Cutoff 0 keeps every cell's top prediction and its confidence, so other cutoffs can be rebuilt offline.
# Builds the augmented h5ads (skipped if already built), then runs one nextflow job.
# Everything goes to cosmos. Run inside tmux; the log goes to logs/.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."
mkdir -p logs
exec > >(tee "logs/test_unlabeled_mm_cutoff${1:-0.25}.log") 2>&1

# The build keeps all unlabeled cells per sample (native ratios). Ratios per study: /space/grp/rschwartz/rschwartz/mask-annotation-overlap/results/qc_labeled_vs_unlabeled/unlabeled_ratio_by_study.tsv (scripts/qc_labeled_vs_unlabeled.py, commit c859220)
cutoff=${1:-0.25}
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

# symlink the cosmos results into the worktree (gitignored)
mkdir -p "$outdir"
ln -sfn "$base/results/census-map-fixes_with_unlabeled_test" results_unlabeled_test

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
    --use_gap "$use_gap" \
    --git_branch census-map-fixes  # reuse the reference caches built on this branch name

find "$outdir" -name '*.unlabeled_unknown_rate.*.tsv.gz'
