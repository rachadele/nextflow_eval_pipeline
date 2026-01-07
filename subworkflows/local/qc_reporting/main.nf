/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC_REPORTING subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Quality control plots and MultiQC report generation
----------------------------------------------------------------------------------------
*/

include { PLOT_QC_COMBINED } from '../../../modules/local/plot_qc_combined/main'
include { MULTIQC          } from '../../../modules/local/multiqc/main'

workflow QC_REPORTING {
    take:
    predicted_meta_channel
    raw_queries_adata

    main:
    // Parse raw queries
    raw_queries_adata_map = raw_queries_adata.map { query_path ->
        def query_name = query_path.getName().split('_raw.h5ad')[0]
        def study_name = query_name.split('_')[0]
        [query_name, study_name, query_path]
    }

    // Parse predicted meta
    predicted_meta_combined = predicted_meta_channel.map { method, query_path, ref_path, predicted_meta ->
        def query_name = query_path.getName().split('.obs.relabel.tsv')[0]
        def study_name = query_name.split('_')[0]
        def ref_name = ref_path.getName().split('\\.')[0]
        [query_name, study_name, method, ref_name, predicted_meta]
    }

    // Combine and group by study/method/ref
    qc_channel_combined = predicted_meta_combined.combine(raw_queries_adata_map, by: [0, 1])

    qc_channel_grouped = qc_channel_combined
        .map { query_name, study_name, method, ref_name, predicted_meta_path, processed_h5ad_path ->
            tuple(
                [study_name, method, ref_name],
                [query_name, predicted_meta_path, processed_h5ad_path]
            )
        }
        .groupTuple(by: 0)
        .map { key, samples ->
            def study_name = key[0]
            def method = key[1]
            def ref_name = key[2]
            [study_name, method, ref_name, samples]
        }

    // Generate QC plots
    PLOT_QC_COMBINED(qc_channel_grouped)

    // Generate MultiQC report
    MULTIQC(PLOT_QC_COMBINED.out.qc_result)

    emit:
    multiqc_html = MULTIQC.out.multiqc_html
}
