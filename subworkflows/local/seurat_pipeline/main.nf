/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SEURAT_PIPELINE subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Seurat label transfer prediction pipeline
----------------------------------------------------------------------------------------
*/

include { QUERY_PROCESS_SEURAT } from '../../../modules/local/query_process_seurat/main'
include { PREDICT_SEURAT       } from '../../../modules/local/predict_seurat/main'
include { CLASSIFY_ALL         } from '../../../modules/local/classify_all/main'

workflow SEURAT_PIPELINE {
    take:
    processed_queries_adata
    ref_paths_seurat
    ref_keys
    ref_region_mapping

    main:
    // Process queries for Seurat
    processed_queries_seurat = QUERY_PROCESS_SEURAT(processed_queries_adata)

    // Combine queries with references
    combos_seurat = processed_queries_seurat.combine(ref_paths_seurat)

    // Run Seurat label transfer
    PREDICT_SEURAT(combos_seurat, ref_keys)

    // Prepare channel for classification
    seurat_scores_channel = PREDICT_SEURAT.out.pred_scores_channel.map { query_path, ref_path, scores_path ->
        ["seurat", query_path, ref_path, scores_path]
    }

    // Classify and compute metrics
    CLASSIFY_ALL(ref_keys, seurat_scores_channel, ref_region_mapping)

    emit:
    f1_scores      = CLASSIFY_ALL.out.f1_score_channel
    predicted_meta = CLASSIFY_ALL.out.predicted_meta_channel
}
