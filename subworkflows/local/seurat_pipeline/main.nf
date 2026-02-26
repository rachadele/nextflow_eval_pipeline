/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SEURAT_PIPELINE subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Seurat label transfer prediction pipeline
----------------------------------------------------------------------------------------
*/

include { QUERY_PROCESS_SEURAT } from "$projectDir/modules/local/query_process_seurat/main"
include { PREDICT_SEURAT       } from "$projectDir/modules/local/predict_seurat/main"
include { CLASSIFY_ALL         } from "$projectDir/modules/local/classify_all/main"

workflow SEURAT_PIPELINE {
    take:
    processed_queries_adata
    ref_paths_seurat
    ref_keys
    ref_region_mapping
    ref_counts

    main:
    // Process queries for Seurat
    processed_queries_seurat = QUERY_PROCESS_SEURAT(processed_queries_adata)

    // Combine queries with references
    combos_seurat = processed_queries_seurat.combine(ref_paths_seurat)

    // Run Seurat label transfer
    PREDICT_SEURAT(combos_seurat, ref_keys)

    // Key ref_counts by ref name for joining
    ref_counts_keyed = ref_counts.map { tsv ->
        [tsv.getName().replace('.ref_counts.tsv', ''), tsv]
    }

    // Prepare channel for classification, joining ref_counts by ref name
    seurat_scores_channel = PREDICT_SEURAT.out.pred_scores_channel
        .map { query_path, ref_path, scores_path ->
            def ref_name = ref_path.toString().split('/').last().replace('.rds', '')
            [ref_name, query_path, ref_path, scores_path]
        }
        .combine(ref_counts_keyed, by: 0)
        .map { ref_name, query_path, ref_path, scores_path, ref_counts_path ->
            ["seurat", query_path, ref_path, scores_path, ref_counts_path]
        }

    // Classify and compute metrics
    CLASSIFY_ALL(ref_keys, seurat_scores_channel, ref_region_mapping)

    emit:
    f1_scores      = CLASSIFY_ALL.out.f1_score_channel
    predicted_meta = CLASSIFY_ALL.out.predicted_meta_channel
}
