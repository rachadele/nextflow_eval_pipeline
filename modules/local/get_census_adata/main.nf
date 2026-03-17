process GET_CENSUS_ADATA {
    label 'process_high'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    storeDir "/cosmos/data/nextflow-eval-pipeline/results/cache/refs/${params.organism}/${params.census_version}/${params.organ}/${params.ref_split}/sub_${params.subsample_ref}"
    publishDir "${params.outdir}/refs/scvi", mode: 'copy', pattern: "**_umap.png"
    publishDir "${params.outdir}/refs/scvi", mode: 'copy', pattern: "refs/*yaml"

    input:
    val ref_collections

    output:
    path "refs/*.h5ad"              , emit: ref_paths_adata
    path "refs/*.ref_counts.tsv"    , emit: ref_counts_adata
    path "**_umap.png"
    path "refs/*yaml"               , emit: ref_region_mapping

    script:
    ref_keys = params.ref_keys.join(' ')
    def author_annotations_arg = params.author_annotations_path ? "--author_annotations_path ${params.author_annotations_path}" : ""
    def original_celltype_columns_arg = params.original_celltype_columns ? "--original_celltype_columns ${params.original_celltype_columns}" : ""
    """
    python $projectDir/bin/get_census_adata.py \\
        --organism ${params.organism} \\
        --census_version ${params.census_version} \\
        --subsample_ref ${params.subsample_ref} \\
        --relabel_path ${params.relabel_r} \\
        --split_column ${params.ref_split} \\
        --ref_collections ${ref_collections} \\
        --ref_keys ${ref_keys} \\
        --organ ${params.organ} \\
        --seed ${params.seed} \\
        ${author_annotations_arg} \\
        ${original_celltype_columns_arg}
    """
}
