#!/bin/bash
# Rebuild per-sample metrics at cutoffs 0 to 0.95 (step 0.05) from the mouse cutoff 0 run. Run inside tmux; log goes to logs/.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."
mkdir -p logs
exec > >(tee logs/rebuild_cutoff_sweep.log) 2>&1
base=/cosmos/data/nextflow-eval-pipeline/results/census-map-fixes_with_unlabeled_test/mus_musculus/ratio_native/ref_500
/home/rschwartz/anaconda3/envs/scanpyenv/bin/python scripts/rebuild_cutoff_metrics.py \
    --validate "$base/cutoff_0.25" --outdir "$base/cutoff_0/rebuilt_cutoffs"
