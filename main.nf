#!/usr/bin/env nextflow

process save_params_to_file {
    publishDir (
        "${params.outdir}",
        mode: "copy"
    )

    output:
    file "params.txt"

    script:
    """
    echo "organism: ${params.organism}" > params.txt
    echo "census_version: ${params.census_version}" >> params.txt
    echo "tree_file: ${params.tree_file}" >> params.txt
    echo "ref_keys: ${params.ref_keys}" >> params.txt
    echo "subsample_ref: ${params.subsample_ref}" >> params.txt
    echo "subsample_query: ${params.subsample_query}" >> params.txt
    echo "relabel_r: ${params.relabel_r}" >> params.txt
    echo "relabel_q: ${params.relabel_q}" >> params.txt
    echo "cutoff: ${params.cutoff}" >> params.txt
    echo "remove_unknown: ${params.remove_unknown}" >> params.txt
    echo "queries_adata: ${params.queries_adata}" >> params.txt
    echo "batch_keys: ${params.batch_keys}" >> params.txt
    echo "ref_split: ${params.ref_split}" >> params.txt
    echo "ref_collections: ${params.ref_collections}" >> params.txt
    echo "integration_method: ${params.integration_method}" >> params.txt
    echo "dims: ${params.dims}" >> params.txt
    echo "max_features: ${params.max_features}" >> params.txt
    echo "k_anchor: ${params.k_anchor}" >> params.txt
    echo "k_score: ${params.k_score}" >> params.txt
    echo "k_weight: ${params.k_weight}" >> params.txt
    echo "outdir: ${params.outdir}" >> params.txt
    """
}


process runSetup {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val organism
    val census_version

    output:
    path "scvi-human-${census_version}/"

    script:
    """
    python $projectDir/bin/setup.py --organism ${organism} --census_version ${census_version}
    """
}

process mapQuery {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val model_path
    tuple val(query_name), val(relabel_q), val(query_file), val(batch_key)
    val ref_keys

    output:
    path "${query_file.getName().toString().replace('.h5ad','_processed.h5ad')}"

script:


"""

python $projectDir/bin/process_query.py \\
                        --model_path ${model_path} \\
                        --relabel_path ${relabel_q} \\
                        --query_path ${query_file} \\
                        --batch_key ${batch_key} \\
                        --subsample_query ${params.subsample_query} \\
                        --ref_keys ${ref_keys} \\
                        ${params.remove_unknown ? '--remove_unknown' : ''}
"""

}

//process processQuerySeurat {
    //conda '/home/rschwartz/anaconda3/envs/r4.3'

    //input:
    //tuple val(query_name), path(relabel_q), path(query_file), val(batch_key)
    //// val subsample_query
    //val ref_keys

    //output:
    //path "${query_file.getName().toString().replace('.rds','_processed.rds')}"


    //script:
    //"""
    //Rscript $projectDir/bin/process_query.R \\
                        //--relabel_path ${relabel_q} \\
                        //--query_path ${query_file} \\
                        //--subsample_query ${params.subsample_query} \\
                        //--ref_keys ${ref_keys}
    //"""
    // 

//}

//process getCensusSeurat {
    //conda '/home/rschwartz/anaconda3/envs/r4.3'

 //input:
    //val organism
    //val census_version
    //val subsample_ref
    //val relabel_r
    //val ref_split
    //val ref_collections

    //output:
    //path "refs/*rds", emit: ref_paths_seurat


    //script:

    //"""
    //# Run the python script to generate the files
    //Rscript $projectDir/bin/get_census_seurat.R \\
        //--organism ${organism} \\
        //--census_version ${census_version} \\
        //--subsample_ref ${subsample_ref} \\
        //--relabel_path ${relabel_r} \\
        //--split_column ${ref_split} \\
        //--ref_collections ${ref_collections}
    //"""
//}


process getCensusAdata {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val organism
    val census_version
    val subsample_ref
    val relabel_r
    val ref_split
    val ref_collections

    output:
    path "refs/*.h5ad", emit: ref_paths_adata

    script:
    """
    # Run the python script to generate the files
    python $projectDir/bin/get_census_adata.py \\
        --organism ${organism} \\
        --census_version ${census_version} \\
        --subsample_ref ${subsample_ref} \\
        --relabel_path ${relabel_r} \\
        --split_column ${ref_split} \\
        --ref_collections ${ref_collections}

    # After running the python script, all .h5ad files will be saved in the refs/ directory inside a work directory
    """
}

process h5adConvertQuery {
    conda '/home/rschwartz/anaconda3/envs/r4.3'

    input:
    path h5ad_file

    output:
    path "${h5ad_file.getName().toString().replace('.h5ad','.rds')}", emit: query_paths_seurat 

    script:
    """
    Rscript $projectDir/bin/h5ad_to_rds.R --h5ad_file ${h5ad_file}
    """
}

process h5adConvertRefs {
    conda '/home/rschwartz/anaconda3/envs/r4.3'

    input:
    path h5ad_file

    output:
    path "${h5ad_file.getName().toString().replace('.h5ad','.rds')}", emit: ref_paths_seurat

    script:
    """
    Rscript $projectDir/bin/h5ad_to_rds.R --h5ad_file ${h5ad_file}
    """
}

process rfPredict {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    // publishDir (
        // path: "${params.outdir}/scvi",
        // mode: "copy"
    // )

    input:
    tuple val(query_path), val(ref_path)
    val ref_keys

    output:
    // path "probs/**"
    tuple path("*obs.relabel.tsv"), val(ref_path), path("probs/*tsv"), emit: probs_channel

    script:
    """
    python $projectDir/bin/rfc_pred.py --query_path ${query_path} --ref_path ${ref_path} --ref_keys ${ref_keys}        
 
    """

}


process predictSeurat {
    conda '/home/rschwartz/anaconda3/envs/r4.3'

    // publishDir (
    //    path "{$params.outdir}/seurat",
    //    mode: "copy"
   // )

    input:
    tuple val(query_path), val(ref_path)
    val ref_keys

    output:
    tuple path("*obs.relabel.tsv"), val(ref_path), path("*prediction_scores_seurat.tsv"), emit: pred_scores_channel 

    script:

    """
    Rscript $projectDir/bin/predict_seurat.R --query_path ${query_path} \\
        --ref_path ${ref_path} \\
        --ref_keys ${ref_keys} \\
        --integration_method ${params.integration_method} \\
        --dims ${params.dims} \\
        --max.features ${params.max_features} \\
        --k.score ${params.k_score} \\
        --k.anchor ${params.k_anchor} \\
        --k.weight ${params.k_weight}
    """

}

process classifyAllAdata {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    publishDir (
        path: "${params.outdir}/scvi",
        mode: "copy"
    )

    input:
    val tree_file
    val ref_keys
    val cutoff
    tuple val(query_path), path(ref_path), path(probs_path)
    val mapping_file

    output:
    path "f1_results/*f1.scores.tsv", emit: f1_score_channel  // Match TSV files in f1_results
    path "roc/**tsv", emit: auc_channel
    path "roc/**png"
    path "confusion/**"
    path "predicted_meta/*tsv"

    script:

    ref_name = ref_path.getName().split('.h5ad')[0]

    """
    python $projectDir/bin/classify_all.py \\
        --tree_file ${tree_file} \\
        --query_path ${query_path} \\
        --ref_name ${ref_name} \\
        --ref_keys ${ref_keys} \\
        --cutoff ${cutoff} \\
        --probs ${probs_path} \\
        --mapping_file ${mapping_file}
    """

}

process classifyAllSeurat {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    publishDir (
        path: "${params.outdir}/seurat",
        mode: "copy"
    )

    input:
    val tree_file
    val ref_keys
    val cutoff
    tuple val(query_path), path(ref_path), path(scores_path)
    val mapping_file

    output:
    path "f1_results/*f1.scores.tsv", emit: f1_score_channel  // Match TSV files in f1_results
    path "roc/**tsv", emit: auc_channel
    path "roc/**png"
    path "confusion/**"
    path "predicted_meta/*tsv"

    script:

    ref_name = ref_path.getName().split('.rds')[0]

    """
    python $projectDir/bin/classify_all.py --tree_file ${tree_file} \\
        --query_path ${query_path} \\
        --ref_name ${ref_name} \\
        --ref_keys ${ref_keys} \\
        --cutoff ${cutoff} \\
        --probs ${scores_path} \\
        --mapping_file ${mapping_file}
    """

}

process plotF1ResultsAdata{
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    publishDir (
        path: "${params.outdir}/scvi",
        mode: "copy"
        )

    input:
    val ref_keys
    val cutoff
    file f1_scores

    output:
    path "f1_plots/*png" // Wildcard to capture all relevant output files
    path "dists/*distribution.png" // Wildcard to capture all relevant output files

    script:
    
    """
    python $projectDir/bin/plot_f1_results.py --ref_keys ${ref_keys} --cutoff ${cutoff} --f1_results ${f1_scores}
 
    """ 
}

process plotF1ResultsSeurat{
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    publishDir (
        path: "${params.outdir}/seurat",
        mode: "copy"
        )

    input:
    val ref_keys
    val cutoff
    file f1_scores

    output:
    path "f1_plots/*png" // Wildcard to capture all relevant output files
    path "dists/*distribution.png" // Wildcard to capture all relevant output files

    script:
    
    """
    python $projectDir/bin/plot_f1_results.py --ref_keys ${ref_keys} --cutoff ${cutoff} --f1_results ${f1_scores}
 
    """ 
}

// Workflow definition
workflow {

    Channel.fromPath(params.queries_adata)
    .set { query_paths_adata }

    Channel.fromPath(params.relabel_q)
    .set { relabel_q_paths }

    // Call the setup process to download the model
    model_path = runSetup(params.organism, params.census_version)
     
    // Get collection names to pull from census
    ref_collections = params.ref_collections.collect { "\"${it}\"" }.join(' ')

    // Get reference data and save to files
    getCensusAdata(params.organism, params.census_version, params.subsample_ref, params.relabel_r, params.ref_split, ref_collections)
    getCensusAdata.out.ref_paths_adata.flatten()
    .set { ref_paths_adata }

    // Convert h5ad files to rds files
     h5adConvertRefs(ref_paths_adata)
     h5adConvertRefs.out.ref_paths_seurat.set { ref_paths_seurat }


    // Get query name from file (including region, eg. Lim_Cingulate)
    relabel_q_paths = relabel_q_paths.map { relabel_q_path -> 
        def relabel_key = relabel_q_path.getName().split('_relabel.tsv')[0]
        [relabel_key, relabel_q_path]
    }

    // Get query names from file (including region)
    query_paths_adata = query_paths_adata.map { query_path -> 
        def query_name = query_path.getName().split('.h5ad')[0]
        [query_name, query_path]
    }

    combined_query_paths_adata = relabel_q_paths
    .join(query_paths_adata)
    .map { query_name, relabel_q_path, query_path ->
        def query_key = query_name.split('_')[0]
        def batch_key = params.batch_keys[query_key]
        [query_key, relabel_q_path, query_path, batch_key]
    }

    //Channel.fromPath(params.queries_seurat)
    //.set { query_paths_seurat }

    //getCensusSeurat(params.organism, params.census_version, params.subsample_ref, params.relabel_r, params.ref_split, ref_collections)
    //getCensusSeurat.out.ref_paths_seurat.flatten()
    //.set { ref_paths_seurat }

    // Map relabeling files to query files
    //query_paths_seurat = query_paths_seurat.map { query_path -> 
        //def query_name = query_path.getName().split('.rds')[0]
        //[query_name, query_path]
    //}

    // Combine the query paths with the relabel paths
    // Map to batch key parameter
    // Returns a tuple with the query prefix (author name only, e.g. Lim), paths and batch key for scvi model


    //combined_query_paths_seurat = relabel_q_paths
    //.join(query_paths_seurat)
    //.map { query_name, relabel_q_path, query_path ->
        //def query_key = query_name.split('_')[0]
        //def batch_key = params.batch_keys[query_key]
        //[query_key, relabel_q_path, query_path, batch_key]
    //}

    // processed_queries_seurat = processQuerySeurat(combined_query_paths_seurat, params.ref_keys.join(' '))

    // Process each query by relabeling, subsampling, and passing through scvi model
    processed_queries_adata = mapQuery(model_path, combined_query_paths_adata, params.ref_keys.join(' ')) 
    // Process each query by relabeling and subsampling
    processed_queries_seurat = h5adConvertQuery(processed_queries_adata)

    // Combine the processed queries with the reference paths
    combos_adata = processed_queries_adata.combine(ref_paths_adata)
    combos_seurat = processed_queries_seurat.combine(ref_paths_seurat)
    
    // Process each query-reference pair
    rfPredict(combos_adata, params.ref_keys.join(' '))
    predictSeurat(combos_seurat, params.ref_keys.join(' '))

    // Collect predictions from each query reference pair
    adata_probs_channel = rfPredict.out.probs_channel 
    seurat_scores_channel = predictSeurat.out.pred_scores_channel


    // Classify all cells based on prediction scores at most granular level
    classifyAllAdata(params.tree_file, params.ref_keys.join(' '), params.cutoff, adata_probs_channel, params.relabel_r)
    f1_scores_adata = classifyAllAdata.out.f1_score_channel

    classifyAllSeurat(params.tree_file, params.ref_keys.join(' '), params.cutoff, seurat_scores_channel, params.relabel_r)
    f1_scores_seurat = classifyAllSeurat.out.f1_score_channel

    // Flatten f1 scores files into a list
    f1_scores_adata 
    .toList()
    .set { f1_scores_adata_files }

    f1_scores_seurat
    .toList()
    .set { f1_scores_seurat_files }

    // Plot f1 score heatmaps using a list of file names from the f1 score channel
    plotF1ResultsAdata(params.ref_keys.join(' '), params.cutoff, f1_scores_adata_files)
    plotF1ResultsSeurat(params.ref_keys.join(' '), params.cutoff, f1_scores_seurat_files)

    save_params_to_file()
}

workflow.onComplete {
    println "Successfully completed"
    println ( workflow.success ? 
    """
    ===============================================================================
    Pipeline execution summary
    -------------------------------------------------------------------------------

    Run as      : ${workflow.commandLine}
    Started at  : ${workflow.start}
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    Config files: ${workflow.configFiles}
    exit status : ${workflow.exitStatus}

    --------------------------------------------------------------------------------
    ================================================================================
    """.stripIndent() : """
    Failed: ${workflow.errorReport}
    exit status : ${workflow.exitStatus}
    """.stripIndent()
    )
}

workflow.onError = {
println "Error: something went wrong, check the pipeline log at '.nextflow.log"
}