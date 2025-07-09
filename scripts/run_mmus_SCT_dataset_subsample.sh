#!/bin/bash

set -e  # Exit on error

# Define parameter values
subsample_ref_values=(50 100 500)
subsample_query=100
ref_split_values=("dataset_id")
cutoff_values=(0 0.05 0.1 0.15 0.2 0.25 0.5 0.75)
normalization_method="SCT"
# not sure if I should handle this now or later
# Loop over parameter combinations
for subsample_ref in "${subsample_ref_values[@]}"; do
    for ref_split in "${ref_split_values[@]}"; do
        for cutoff in "${cutoff_values[@]}"; do
                echo "Running: subsample_ref=$subsample_ref, ref_split=$ref_split, cutoff=$cutoff, normalization_method=$normalization_method"
                nextflow main.nf -params-file /space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/params.mm.json \
                    --subsample_query "$subsample_query" \
                    --subsample_ref "$subsample_ref" \
                    --ref_split "$ref_split" \
                    -profile conda \
                    --cutoff "$cutoff" \
                    --subset_type sample \
                    -work-dir mmus \
                    --batch_correct true \
                    -resume \
                    --remove_unknown \
                    --normalization_method "$normalization_method" \
                    --outdir_prefix "mus_musculus/subsample_query_100/old_refs" \
                    -process.executor "slurm"
            done
        done
    done
