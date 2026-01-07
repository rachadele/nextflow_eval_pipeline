# Copilot Instructions for nextflow_eval_pipeline

## Project Overview
- This project benchmarks automated cell-type annotation (label transfer) for single-cell RNA-seq (scRNA-seq) data, focusing on human and mouse neocortex datasets.
- Two main annotation strategies are evaluated: Random Forest on SCVI embeddings (Python) and Seurat label transfer (R).
- The pipeline is orchestrated with Nextflow (`main.nf`), integrating Python and R scripts for data processing, model training, and evaluation.

## Key Components
- **main.nf**: Nextflow pipeline orchestrating all steps. Entry point for all workflows.
- **bin/**: Contains all main scripts:
  - Python: `classify_all.py`, `predict_scvi.py`, `get_census_adata.py`, `process_query.py`, etc.
  - R: `predict_seurat.R`, `seurat_preprocessing.R`, `ref_preprocessing.R`, etc.
- **meta/**: Reference and mapping files (e.g., `census_map_human.tsv`, `cell_type_markers.tsv`).
- **params.*.json**: Example parameter sets for different runs.
- **nextflow.config**: All pipeline parameters, including reference/test data, batch keys, and method settings.

## Developer Workflows
- **Run the pipeline**: `nextflow run main.nf -profile local` (adjust parameters in `nextflow.config` or via CLI).
- **Add new scripts**: Place in `bin/` and update `main.nf` to add new processes.
- **Reference data**: Update or add mapping/marker files in `meta/`.
- **Environment**: Python and R environments are managed via conda (see `main.nf` for env paths).

## Project-Specific Conventions
- **Label keys**: Always use `subclass`, `class`, `family`, `global` as main cell type levels (see `params.ref_keys`).
- **Reference/test splits**: Controlled by `params.ref_collections`, `params.ref_split`, and `params.subsample_*` in `nextflow.config`.
- **Batch correction**: Use `params.batch_correct` for Seurat integration.
- **Output structure**: Results are organized by census version, organism, and parameterization (see `params.outdir`).
- **Random seeds**: Set in both Python and R scripts for reproducibility (default: 42).

## Integration & Data Flow
- Data flows from reference/test data (in `2024-07-01/`, `mmus/`, `hsap/`, etc.) through preprocessing, model training, prediction, and evaluation.
- Cross-language integration: Nextflow manages handoff between Python and R scripts, passing files as process outputs/inputs.
- All scripts expect explicit CLI arguments; see script headers for details.

## Examples
- To add a new annotation method, create a script in `bin/`, add a Nextflow process, and update `nextflow.config` as needed.
- To benchmark a new dataset, add it to the appropriate data folder and update `params.*.json` or `nextflow.config`.

## References
- See `README.md` for scientific background and references.
- See `nextflow.config` for all tunable parameters and conventions.

---
For questions about workflow structure, see `main.nf` and `README.md`. For parameter conventions, see `nextflow.config`.
