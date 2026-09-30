#!/bin/bash
# Mouse cutoff sweep on queries that include author-unlabeled cells
# (built by scripts/build_unlabeled_queries.py). Same settings as
# run_mmus_sweep_census_fixes.sh, but reads the augmented h5ads and writes to a
# separate outdir so the census-fix sweep results are not overwritten.
# Run from the worktree root. USE_GAP=true sweeps the confidence-gap cutoff instead.
set -e
queries=/space/grp/rschwartz/rschwartz/get_gemma_data.nf/study_names_mouse.txt_author_true_process_samples_true_with_unlabeled/h5ad/**h5ad
results=/cosmos/data/nextflow-eval-pipeline/results/census-map-fixes_with_unlabeled
use_gap=${USE_GAP:-false}
subsample_ref_values=(${SUBSAMPLE_REFS:-500 100 50})
subsample_query=100
ref_split_values=("dataset_id")
cutoff_values=(${CUTOFFS:-0 0.05 0.1 0.15 0.2 0.25 0.5 0.75})
normalization_method="SCT"

for subsample_ref in "${subsample_ref_values[@]}"; do
    for ref_split in "${ref_split_values[@]}"; do
        for cutoff in "${cutoff_values[@]}"; do
            outdir="$results/2024-07-01/mus_musculus/$subsample_query/$ref_split/$normalization_method/gap_$use_gap/ref_$subsample_ref/cutoff_$cutoff"
            echo "=== Running: subsample_ref=$subsample_ref, ref_split=$ref_split, cutoff=$cutoff, use_gap=$use_gap -> $outdir ==="
            nextflow main.nf -params-file params.mm.json \
                --queries_adata "$queries" \
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
        done
    done
done
echo "=== SWEEP COMPLETE ==="
