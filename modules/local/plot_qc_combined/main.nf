process PLOT_QC_COMBINED {
    label 'process_medium'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/${method}/${study_name}/${ref_name}/qc_plots", mode: 'copy'

    input:
    tuple val(study_name), val(method), val(ref_name), val(samples)

    output:
    tuple val(study_name), val(method), val(ref_name), path("${study_name}/"), emit: qc_result
    path "**_mqc.png"
    path "**.tsv"
    path "**.txt"

    script:
    predicted_meta_list = samples.collect { it[1] }.join(" ")
    query_path_list = samples.collect { it[2] }.join(" ")
    """
    python $projectDir/bin/plot_QC_combined.py \\
        --query_paths ${query_path_list} \\
        --predicted_meta_files ${predicted_meta_list} \\
        --organism ${params.organism} \\
        --markers_file ${params.markers_file} \\
        --nmads 5 \\
        --gene_mapping ${params.gene_mapping} \\
        --ref_keys ${params.ref_keys.join(' ')} \\
        --mapping_file ${params.relabel_r} \\
        --study_name ${study_name}
    """
}
