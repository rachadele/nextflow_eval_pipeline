/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCVI_PIPELINE subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCVI-based Random Forest prediction pipeline
----------------------------------------------------------------------------------------
*/

include { RF_PREDICT   } from '../../../modules/local/rf_predict/main'
include { CLASSIFY_ALL } from '../../../modules/local/classify_all/main'

workflow SCVI_PIPELINE {
    take:
    processed_queries_adata
    ref_paths_adata
    ref_keys
    ref_region_mapping

    main:
    // Combine queries with references
    combos_adata = processed_queries_adata.combine(ref_paths_adata)

    // Run RF prediction on SCVI embeddings
    RF_PREDICT(combos_adata, ref_keys)

    // Prepare channel for classification
    adata_probs_channel = RF_PREDICT.out.probs_channel.map { query_path, ref_path, probs_path ->
        ["scvi", query_path, ref_path, probs_path]
    }

    // Classify and compute metrics
    CLASSIFY_ALL(ref_keys, adata_probs_channel, ref_region_mapping)

    emit:
    f1_scores      = CLASSIFY_ALL.out.f1_score_channel
    predicted_meta = CLASSIFY_ALL.out.predicted_meta_channel
}
