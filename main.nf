#!/usr/bin/env nextflow

process save_params_to_file {
    publishDir (
        "${params.outdir}",
        mode: "copy"
    )

    output:
    file "params.yaml"

    script:
    """
    cat <<EOF > params.yaml
    organism: ${params.organism}
    census_version: ${params.census_version}
    ref_keys: ${params.ref_keys}
    subsample_ref: ${params.subsample_ref}
    subsample_query: ${params.subsample_query}
    relabel_r: ${params.relabel_r}
    relabel_q: ${params.relabel_q}
    cutoff: ${params.cutoff}
    remove_unknown: ${params.remove_unknown}
    queries_adata: ${params.queries_adata}
    batch_keys: ${params.batch_keys}
    ref_split: ${params.ref_split}
    ref_collections: ${params.ref_collections}
    integration_method: ${params.integration_method}
    dims: ${params.dims}
    max_features: ${params.max_features}
    k_anchor: ${params.k_anchor}
    k_score: ${params.k_score}
    k_weight: ${params.k_weight}
    outdir: ${params.outdir}
    normalization_method: ${params.normalization_method}
    subset_type: ${params.subset_type}
    batch_correct: ${params.batch_correct}
    EOF
    """
}



process runSetup {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val organism
    val census_version

    output:
    path "scvi-${params.organism}-${census_version}/"

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
                            ${params.subsample_query != null ? "--subsample_query ${params.subsample_query}" : ""} \
                            --ref_keys ${ref_keys} \\
                            --seed ${params.seed} \\
                            ${params.remove_unknown ? '--remove_unknown' : ''}
    """

}

process getCensusAdata {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    publishDir "${params.outdir}", pattern: "**obs.tsv", mode: "copy"

    input:
    val organism
    val census_version
    val subsample_ref
    val relabel_r
    val ref_split
    val ref_collections

    output:
    path "refs/*.h5ad", emit: ref_paths_adata
    path "**obs.tsv"
    path "refs/*yaml", emit: ref_region_mapping

    script:
    ref_keys = params.ref_keys.join(' ')
    """
    # Run the python script to generate the files
    python $projectDir/bin/get_census_adata.py \\
        --organism ${organism} \\
        --census_version ${census_version} \\
        --subsample_ref ${subsample_ref} \\
        --relabel_path ${relabel_r} \\
        --split_column ${ref_split} \\
        --ref_collections ${ref_collections} \\
        --ref_keys ${ref_keys} \\
        --seed ${params.seed}

    # After running the python script, all .h5ad files will be saved in the refs/ directory inside a work directory
    """
}

process queryProcessSeurat {
    conda '/home/rschwartz/anaconda3/envs/r4.3'

    input:
    path h5ad_file

    output:
    path "${h5ad_file.getName().toString().replace('.h5ad','.rds')}", emit: query_paths_seurat 

    script:
    """
    Rscript $projectDir/bin/seurat_preprocessing.R --h5ad_file ${h5ad_file} \\
            --normalization_method ${params.normalization_method} \\
            --dims ${params.dims} \\
            --nfeatures ${params.nfeatures}

    """
}

process refProcessSeurat {
    conda '/home/rschwartz/anaconda3/envs/r4.3'

    input:
    path h5ad_file

    output:
    path "${h5ad_file.getName().toString().replace('.h5ad','.rds')}", emit: ref_paths_seurat

    script:
    """
    Rscript $projectDir/bin/ref_preprocessing.R --h5ad_file ${h5ad_file} \\
            --normalization_method ${params.normalization_method} \\
            --dims ${params.dims} \\
            --nfeatures ${params.nfeatures} \\
            ${params.batch_correct ? '--batch_correct' : ''} 
            
    """
}

process rfPredict {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    tuple val(query_path), val(ref_path)
    val ref_keys

    output:
    tuple path("*obs.relabel.tsv"), val(ref_path), path("probs/*tsv"), emit: probs_channel

    script:
    """
    python $projectDir/bin/predict_scvi.py --query_path ${query_path} --ref_path ${ref_path} --ref_keys ${ref_keys}        
 
    """

}


process predictSeurat {
    conda '/home/rschwartz/anaconda3/envs/r4.3'


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
        --k.weight ${params.k_weight} \\
        --normalization_method ${params.normalization_method}
    """

}

process classifyAllAdata {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    publishDir (
        path: "${params.outdir}/scvi",
        mode: "copy"
    )

    input:
    val ref_keys
    val cutoff
    tuple val(query_path), path(ref_path), path(probs_path)
    val mapping_file
    val ref_region_mapping

    output:
    path "f1_results/*f1.scores.tsv", emit: f1_score_channel  // Match TSV files in f1_results
   // path "roc/**tsv", emit: auc_channel
    //path "roc/**png"
    path "confusion/**"
    path "predicted_meta/*tsv"

    script:

    ref_name = ref_path.getName().split('.h5ad')[0]

    """
    python $projectDir/bin/classify_all.py \\
        --query_path ${query_path} \\
        --ref_name ${ref_name} \\
        --ref_keys ${ref_keys} \\
        --cutoff ${cutoff} \\
        --probs ${probs_path} \\
        --mapping_file ${mapping_file} \\
        --ref_region_mapping ${ref_region_mapping}
    """

}

process classifyAllSeurat {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    publishDir (
        path: "${params.outdir}/seurat",
        mode: "copy"
    )

    input:
    val ref_keys
    val cutoff
    tuple val(query_path), path(ref_path), path(scores_path) 
    val mapping_file
    val ref_region_mapping

    output:
    path "f1_results/*f1.scores.tsv", emit: f1_score_channel  // Match TSV files in f1_results
    //path "roc/**tsv", emit: auc_channel
   // path "roc/**png"
    path "confusion/**"
    path "predicted_meta/*tsv"

    script:

    ref_name = ref_path.getName().split('.rds')[0]

    """
    python $projectDir/bin/classify_all.py  \\
        --query_path ${query_path} \\
        --ref_name ${ref_name} \\
        --ref_keys ${ref_keys} \\
        --cutoff ${cutoff} \\
        --probs ${scores_path} \\
        --mapping_file ${mapping_file} \\
        --ref_region_mapping ${ref_region_mapping}
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
    getCensusAdata.out.ref_region_mapping.set { ref_region_mapping }
    // Convert h5ad files to rds files
     refProcessSeurat(ref_paths_adata)
     refProcessSeurat.out.ref_paths_seurat.set { ref_paths_seurat }



    // Get query name from file (including region, eg. Lim_Cingulate)
    relabel_q_paths = relabel_q_paths.map { relabel_q_path -> 
        def query_key = relabel_q_path.getName().split('_relabel.tsv')[0]
        [query_key, relabel_q_path]
    }

    // Get query names from file (including region)
    query_paths_adata = query_paths_adata.map { query_path -> 
        def query_name = query_path.getName().split('.h5ad')[0]
        def query_key = query_name.split('_')[0]
        [query_key, query_name, query_path]
    }
 
    combined_query_paths_adata = query_paths_adata
        .combine(relabel_q_paths, by: 0) // Match query_key
        .map { query_key, query_name, query_path, relabel_q_path ->
            def batch_key = params.batch_keys.get(query_key, "sample_id")  // Use default if not found
            [query_name, relabel_q_path, query_path, batch_key]
        }
    // Process each query by relabeling, subsampling, and passing through scvi model
    processed_queries_adata = mapQuery(model_path, combined_query_paths_adata, params.ref_keys.join(' ')) 
    // Process each query by relabeling and subsampling
    processed_queries_seurat = queryProcessSeurat(processed_queries_adata)

    // Combine the processed queries with the reference paths
    combos_adata = processed_queries_adata.combine(ref_paths_adata)
    combos_seurat = processed_queries_seurat.combine(ref_paths_seurat)
    // Process each query-reference pair
    rfPredict(combos_adata, params.ref_keys.join(' '))
    predictSeurat(combos_seurat, params.ref_keys.join(' '))

    // Collect predictions from each query reference pair
    adata_probs_channel = rfPredict.out.probs_channel 
    seurat_scores_channel = predictSeurat.out.pred_scores_channel


    adata_probs_channel = adata_probs_channel.map { query_path, ref_path, probs_path ->
        query_name = query_path.getName().split('.obs.relabel.tsv')[0]
        //query_region = params.query_region[query_name]
       // query_sex = params.query_sex[query_name]
       // query_disease = params.query_disease[query_name]
        [query_path, ref_path, probs_path]
    }
    seurat_scores_channel = seurat_scores_channel.map { query_path, ref_path, scores_path ->
        query_name = query_path.getName().split('.obs.relabel.tsv')[0]
        //query_region = params.query_region[query_name]
       // query_sex = params.query_sex[query_name]
       // query_disease = params.query_disease[query_name]
        [query_path, ref_path, scores_path]
    }
    // Classify all cells based on prediction scores at most granular level
    classifyAllAdata(params.ref_keys.join(' '), params.cutoff, adata_probs_channel, params.relabel_r, ref_region_mapping)
    f1_scores_adata = classifyAllAdata.out.f1_score_channel

    classifyAllSeurat(params.ref_keys.join(' '), params.cutoff, seurat_scores_channel, params.relabel_r, ref_region_mapping)
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
