process SAVE_PARAMS {
    label 'process_low'

    output:
    path "params.yaml"

    script:
    """
    cat <<EOF > params.yaml
    organism: ${params.organism}
    census_version: ${params.census_version}
    ref_keys: ${params.ref_keys}
    subsample_ref: ${params.subsample_ref}
    subsample_query: ${params.subsample_query}
    relabel_r: ${params.relabel_r}
    relabel_q: ${params.relabel_q}
    cutoff: ${params.cutoff}
    remove_unknown: ${params.remove_unknown}
    queries_adata: ${params.queries_adata}
    batch_keys: ${params.batch_keys}
    ref_split: ${params.ref_split}
    ref_collections: ${params.ref_collections}
    integration_method: ${params.integration_method}
    dims: ${params.dims}
    max_features: ${params.max_features}
    k_anchor: ${params.k_anchor}
    k_score: ${params.k_score}
    k_weight: ${params.k_weight}
    outdir: ${params.outdir}
    normalization_method: ${params.normalization_method}
    subset_type: ${params.subset_type}
    batch_correct: ${params.batch_correct}
    nmads: ${params.nmads}
    use_gap: ${params.use_gap}
    EOF
    """
}
