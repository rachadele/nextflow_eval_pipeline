process RF_PREDICT {
    label 'process_medium'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    tuple val(query_path), val(ref_path)
    val ref_keys

    output:
    tuple path("*obs.relabel.tsv"), val(ref_path), path("probs/*tsv"), emit: probs_channel

    script:
    ref_name = ref_path.getName().split('\\.h5ad')[0]
    query_name = query_path.getName().split('\\_processed.h5ad')[0]
    study_name = query_name.split('_')[0]
    """
    python $projectDir/bin/predict_scvi.py \\
        --query_path ${query_path} \\
        --ref_path ${ref_path} \\
        --ref_keys ${ref_keys}
    """
}
