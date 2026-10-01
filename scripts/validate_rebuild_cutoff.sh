#!/bin/bash
# Check that metrics rebuilt from the cutoff 0 run reproduce the real cutoff 0.25 run. Run inside tmux; log goes to logs/.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."
mkdir -p logs
exec > >(tee logs/validate_rebuild_cutoff.log) 2>&1
base=/cosmos/data/nextflow-eval-pipeline/results/census-map-fixes_with_unlabeled_test/mus_musculus/ratio_native/ref_500
/home/rschwartz/anaconda3/envs/scanpyenv/bin/python scripts/rebuild_cutoff_metrics.py --cutoffs 0.25 \
    --validate "$base/cutoff_0.25" --outdir "$base/cutoff_0/rebuilt_cutoffs_validation" "$@"
