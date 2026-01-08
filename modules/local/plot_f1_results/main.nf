process PLOT_F1_RESULTS_SCVI {
    label 'process_low'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/scvi", mode: 'copy'

    input:
    val ref_keys
    val cutoff
    path f1_scores

    output:
    path "dists/*distribution.png"

    script:
    """
    python $projectDir/bin/plot_results_summary.py \\
        --ref_keys ${ref_keys} \\
        --cutoff ${cutoff} \\
        --f1_results ${f1_scores}
    """
}

process PLOT_F1_RESULTS_SEURAT {
    label 'process_low'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/seurat", mode: 'copy'

    input:
    val ref_keys
    val cutoff
    path f1_scores

    output:
    path "dists/*distribution.png"

    script:
    """
    python $projectDir/bin/plot_results_summary.py \\
        --ref_keys ${ref_keys} \\
        --cutoff ${cutoff} \\
        --f1_results ${f1_scores}
    """
}
