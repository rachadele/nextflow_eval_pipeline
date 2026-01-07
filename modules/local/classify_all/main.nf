process CLASSIFY_ALL {
    label 'process_medium'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val ref_keys
    tuple val(method), path(query_path), path(ref_path), path(probs_path)
    val ref_region_mapping

    output:
    tuple val(method), path("label_transfer_metrics/*summary.scores.tsv"), emit: f1_score_channel
    path "confusion/**"
    tuple val(method), path("${query_path}"), path("${ref_path}"), path("predicted_meta/**tsv"), emit: predicted_meta_channel

    script:
    def ref_name = ref_path.getName().split('\\.')[0]
    def query_name = query_path.getName().split('\\.obs.relabel.tsv')[0]
    def study_name = query_name.split('_')[0]
    """
    python $projectDir/bin/classify_all.py \\
        --query_path ${query_path} \\
        --ref_name ${ref_name} \\
        --ref_keys ${ref_keys} \\
        --cutoff ${params.cutoff} \\
        --probs ${probs_path} \\
        --mapping_file ${params.relabel_r} \\
        --ref_region_mapping ${ref_region_mapping} \\
        --study_name ${study_name} \\
        --method ${method} \\
        ${params.use_gap ? '--use_gap' : ''}
    """
}
