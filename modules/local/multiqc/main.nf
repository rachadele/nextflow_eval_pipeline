process MULTIQC {
    label 'process_low'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/${method}/${study_name}/${ref_name}", mode: 'copy', pattern: "**multiqc_report.html"

    input:
    tuple val(study_name), val(method), val(ref_name), path(qc_dir)

    output:
    tuple val(study_name), val(method), val(ref_name), path("**multiqc_report.html"), emit: multiqc_html

    script:
    """
    cp ${params.multiqc_config} new_config.yaml
    echo 'title: "${params.census_version} ${study_name} ${method} ${ref_name} ${params.cutoff}"' >> new_config.yaml
    multiqc ${qc_dir} --config new_config.yaml
    """
}
