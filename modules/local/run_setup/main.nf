process RUN_SETUP {
    label 'process_low'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    storeDir "/cosmos/data/nextflow-eval-pipeline/results/cache/scvi_model/${params.organism}/${census_version}"

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
