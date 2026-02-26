/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCVI_PIPELINE subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCVI-based Random Forest prediction pipeline
----------------------------------------------------------------------------------------
*/

include { RF_PREDICT   } from "$projectDir/modules/local/rf_predict/main"
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

    // Run RF prediction on SCVI embeddings
    RF_PREDICT(combos_adata, ref_keys)

    // Key ref_counts by ref name for joining
    ref_counts_keyed = ref_counts.map { tsv ->
        [tsv.getName().replace('.ref_counts.tsv', ''), tsv]
    }

    // Prepare channel for classification, joining ref_counts by ref name
    adata_probs_channel = RF_PREDICT.out.probs_channel
        .map { query_path, ref_path, probs_path ->
            def ref_name = ref_path.toString().split('/').last().replace('.h5ad', '')
            [ref_name, query_path, ref_path, probs_path]
        }
        .combine(ref_counts_keyed, by: 0)
        .map { ref_name, query_path, ref_path, probs_path, ref_counts_path ->
            ["scvi", query_path, ref_path, probs_path, ref_counts_path]
        }

    // Classify and compute metrics
    CLASSIFY_ALL(ref_keys, adata_probs_channel, ref_region_mapping)

    emit:
    f1_scores      = CLASSIFY_ALL.out.f1_score_channel
    predicted_meta = CLASSIFY_ALL.out.predicted_meta_channel
}
