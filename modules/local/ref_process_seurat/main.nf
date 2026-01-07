process REF_PROCESS_SEURAT {
    label 'process_high'
    conda '/home/rschwartz/anaconda3/envs/r4.3'

    input:
    path h5ad_file

    output:
    path "${h5ad_file.getName().toString().replace('.h5ad','.rds')}", emit: ref_paths_seurat

    script:
    """
    Rscript $projectDir/bin/ref_preprocessing.R \\
        --h5ad_file ${h5ad_file} \\
        --normalization_method ${params.normalization_method} \\
        --dims ${params.dims} \\
        --nfeatures ${params.nfeatures} \\
        --k.score ${params.k_score} \\
        --k.anchor ${params.k_anchor} \\
        --k.weight ${params.k_weight} \\
        ${params.batch_correct ? '--batch_correct' : ''}
    """
}
