#!/bin/bash

set -e  # Exit on error

# Define parameter values
subsample_ref_values=(500)
census_versions=("2025-01-30")
#subsample_query=100
ref_split_values=("dataset_id")
cutoff_values=(0)
normalization_method="SCT"
# Loop over parameter combinations
for subsample_ref in "${subsample_ref_values[@]}"; do
    for ref_split in "${ref_split_values[@]}"; do
        for cutoff in "${cutoff_values[@]}"; do
			for census_version in "${census_versions[@]}"; do
                echo "Running: subsample_ref=$subsample_ref, ref_split=$ref_split, cutoff=$cutoff, normalization_method=$normalization_method"
                nextflow main.nf -params-file /space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/params.mm.json \
                    --subsample_ref "$subsample_ref" \
                    --ref_split "$ref_split" \
                    -profile conda \
                    --cutoff "$cutoff" \
                    --subset_type sample \
                    -work-dir mmus_minimal \
                    --batch_correct true \
                    --normalization_method "$normalization_method" \
                    -process.executor slurm \
					--census_version "$census_version" \
                    -resume \
                    --remove_unknown false \
					--outdir_prefix "$census_version/mus_musculus/minimal/keep_unknown" \
                    --use_gap true
				    #--subsample_query "$subsample_query" \
                done
            done
        done
    done