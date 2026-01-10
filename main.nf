#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Cell Type Annotation Benchmarking Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Benchmarks SCVI-based Random Forest vs Seurat label transfer methods
    for automated cell type annotation of scRNA-seq data using CellXGene Census
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAVE_PARAMS        } from "$projectDir/modules/local/save_params/main"
include { MAP_QUERY          } from "$projectDir/modules/local/map_query/main"
include { PREPARE_REFERENCES } from "$projectDir/subworkflows/local/prepare_references/main"
include { SCVI_PIPELINE      } from "$projectDir/subworkflows/local/scvi_pipeline/main"
include { SEURAT_PIPELINE    } from "$projectDir/subworkflows/local/seurat_pipeline/main"
include { QC_REPORTING       } from "$projectDir/subworkflows/local/qc_reporting/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // Set up input channels
    query_paths_adata = Channel.fromPath(params.queries_adata)
    relabel_q_paths = Channel.fromPath(params.relabel_q)

    // Format collection names for Census query
    ref_collections = params.ref_collections.collect { "\"${it}\"" }.join(' ')

    // Prepare reference data (SCVI model + Census references)
    PREPARE_REFERENCES(
        params.organism,
        params.census_version,
        ref_collections
    )

    // Map relabel paths to query keys
    relabel_q_paths = relabel_q_paths.map { relabel_q_path ->
        def query_key = relabel_q_path.getName().split('_relabel.tsv')[0]
        [query_key, relabel_q_path]
    }

    // Map query paths to query names and keys
    query_paths_adata = query_paths_adata.map { query_path ->
        def query_name = query_path.getName().split('.h5ad')[0]
        def query_key = query_name.split('_')[0]
        [query_key, query_name, query_path]
    }

    // Combine queries with relabel files and batch keys
    combined_query_paths_adata = query_paths_adata
        .combine(relabel_q_paths, by: 0)
        .map { query_key, query_name, query_path, relabel_q_path ->
            def batch_key = params.batch_keys.get(query_key, "sample_id")
            [query_name, relabel_q_path, query_path, batch_key]
        }

    // Process queries through SCVI model
    MAP_QUERY(
        PREPARE_REFERENCES.out.model_path,
        combined_query_paths_adata,
        params.ref_keys.join(' ')
    )

    processed_queries_adata = MAP_QUERY.out.processed_query_adata
    raw_queries_adata = MAP_QUERY.out.raw_query_adata

    // Run SCVI pipeline
    SCVI_PIPELINE(
        processed_queries_adata,
        PREPARE_REFERENCES.out.ref_paths_adata,
        params.ref_keys.join(' '),
        PREPARE_REFERENCES.out.ref_region_mapping
    )

    // Run Seurat pipeline
    SEURAT_PIPELINE(
        processed_queries_adata,
        PREPARE_REFERENCES.out.ref_paths_seurat,
        params.ref_keys.join(' '),
        PREPARE_REFERENCES.out.ref_region_mapping
    )

    // Combine predicted metadata from both methods
    predicted_meta_combined = SCVI_PIPELINE.out.predicted_meta.concat(
        SEURAT_PIPELINE.out.predicted_meta
    )

    // Generate QC reports
    QC_REPORTING(
        predicted_meta_combined,
        raw_queries_adata
    )

    // Save parameters
    SAVE_PARAMS()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    println "Successfully completed"
    println(workflow.success ?
        """
        ================================================================================
        Pipeline execution summary
        --------------------------------------------------------------------------------

        Run as      : ${workflow.commandLine}
        Started at  : ${workflow.start}
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        Config files: ${workflow.configFiles}
        Exit status : ${workflow.exitStatus}
        Output dir  : ${params.outdir}

        --------------------------------------------------------------------------------
        ================================================================================
        """.stripIndent() :
        """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """.stripIndent()
    )
}

workflow.onError = {
    println "Error: something went wrong, check the pipeline log at '.nextflow.log'"
}
