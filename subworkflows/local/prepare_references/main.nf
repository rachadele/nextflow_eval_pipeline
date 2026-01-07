/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_REFERENCES subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downloads and prepares reference data from CellXGene Census
----------------------------------------------------------------------------------------
*/

include { RUN_SETUP           } from '../../../modules/local/run_setup/main'
include { GET_CENSUS_ADATA    } from '../../../modules/local/get_census_adata/main'
include { REF_PROCESS_SEURAT  } from '../../../modules/local/ref_process_seurat/main'

workflow PREPARE_REFERENCES {
    take:
    organism
    census_version
    ref_collections

    main:
    // Download SCVI model
    model_path = RUN_SETUP(organism, census_version)

    // Get reference data from Census
    GET_CENSUS_ADATA(ref_collections)

    ref_paths_adata = GET_CENSUS_ADATA.out.ref_paths_adata.flatten()
    ref_region_mapping = GET_CENSUS_ADATA.out.ref_region_mapping

    // Process references for Seurat
    REF_PROCESS_SEURAT(ref_paths_adata)
    ref_paths_seurat = REF_PROCESS_SEURAT.out.ref_paths_seurat

    emit:
    model_path         = model_path
    ref_paths_adata    = ref_paths_adata
    ref_paths_seurat   = ref_paths_seurat
    ref_region_mapping = ref_region_mapping
}
