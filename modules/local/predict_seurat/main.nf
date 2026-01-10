process PREDICT_SEURAT {
    label 'process_high'
    conda '/home/rschwartz/anaconda3/envs/r4.3'

    input:
    tuple val(query_path), val(ref_path)
    val ref_keys

    output:
    tuple path("*obs.relabel.tsv"), val(ref_path), path("*prediction_scores_seurat.tsv"), emit: pred_scores_channel

    script:
    """
    Rscript $projectDir/bin/predict_seurat.R \\
        --query_path ${query_path} \\
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
