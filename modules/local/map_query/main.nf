process MAP_QUERY {
    label 'process_high'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val model_path
    tuple val(query_name), val(relabel_q), val(query_file), val(batch_key)
    val ref_keys

    output:
    path "${query_file.getName().toString().replace('.h5ad','_processed.h5ad')}", emit: processed_query_adata
    path "${query_file.getName().toString().replace('.h5ad','_raw.h5ad')}"      , emit: raw_query_adata

    script:
    """
    python $projectDir/bin/process_query.py \\
        --model_path ${model_path} \\
        --relabel_path ${relabel_q} \\
        --query_path ${query_file} \\
        --batch_key ${batch_key} \\
        ${params.subsample_query != null ? "--subsample_query ${params.subsample_query}" : ""} \\
        --ref_keys ${ref_keys} \\
        --seed ${params.seed} \\
        --nmads ${params.nmads} \\
        --gene_mapping ${params.gene_mapping} \\
        --mapping_file ${params.relabel_r} \\
        ${params.remove_unknown ? '--remove_unknown' : ''}
    """
}
