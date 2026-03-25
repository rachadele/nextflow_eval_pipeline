/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCVI_PIPELINE subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCVI-based prediction pipeline (Random Forest and kNN classifiers)
----------------------------------------------------------------------------------------
*/

include { SCVI_PREDICT } from "$projectDir/modules/local/scvi_predict/main"
include { CLASSIFY_ALL } from "$projectDir/modules/local/classify_all/main"

workflow SCVI_PIPELINE {
    take:
    processed_queries_adata
    ref_paths_adata
    ref_keys
    ref_region_mapping
    ref_counts

    main:
    // Combine queries with references
    combos_adata = processed_queries_adata.combine(ref_paths_adata)

    // Run both RF and kNN prediction in a single process per query-ref pair
    SCVI_PREDICT(combos_adata, ref_keys)

    // Key ref_counts by ref name for joining
    ref_counts_keyed = ref_counts.map { tsv ->
        [tsv.getName().replace('.ref_counts.tsv', ''), tsv]
    }

    // Split the two prob files into separate channel entries, one per classifier
    adata_probs_channel = SCVI_PREDICT.out.probs_channel
        .map { query_path, ref_path, rf_probs, knn_probs ->
            def ref_name = ref_path.toString().split('/').last().replace('.h5ad', '')
            [ref_name, query_path, ref_path, rf_probs, knn_probs]
        }
        .combine(ref_counts_keyed, by: 0)
        .flatMap { ref_name, query_path, ref_path, rf_probs, knn_probs, ref_counts_path ->
            [
                ["scvi_rf",  query_path, ref_path, rf_probs,  ref_counts_path],
                ["scvi_knn", query_path, ref_path, knn_probs, ref_counts_path]
            ]
        }

    // Classify and compute metrics
    CLASSIFY_ALL(ref_keys, adata_probs_channel, ref_region_mapping)

    emit:
    f1_scores      = CLASSIFY_ALL.out.f1_score_channel
    predicted_meta = CLASSIFY_ALL.out.predicted_meta_channel
}
