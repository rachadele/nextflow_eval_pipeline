# nf-core Reorganization Summary

This document summarizes the changes made to reorganize the pipeline according to nf-core best practices.

---

## Overview

The monolithic `main.nf` (557 lines) was refactored into a modular DSL2 structure with separate modules, subworkflows, and configuration files.

---

## New Directory Structure

```
nextflow_eval_pipeline/
├── main.nf                     # Reduced to ~150 lines (orchestrator only)
├── nextflow.config             # Imports config files, defines profiles
├── conf/
│   ├── base.config             # Process defaults, resource labels
│   ├── params.config           # Default parameter values
│   ├── modules.config          # Per-module publishDir settings
│   ├── test_hsap.config        # Human test parameters
│   └── test_mmus.config        # Mouse test parameters
├── modules/local/
│   ├── save_params/main.nf
│   ├── run_setup/main.nf
│   ├── get_census_adata/main.nf
│   ├── map_query/main.nf
│   ├── query_process_seurat/main.nf
│   ├── ref_process_seurat/main.nf
│   ├── rf_predict/main.nf
│   ├── predict_seurat/main.nf
│   ├── classify_all/main.nf
│   ├── plot_f1_results/main.nf
│   ├── plot_qc_combined/main.nf
│   └── multiqc/main.nf
├── subworkflows/local/
│   ├── prepare_references/main.nf
│   ├── scvi_pipeline/main.nf
│   ├── seurat_pipeline/main.nf
│   └── qc_reporting/main.nf
├── tests/
│   ├── run_test_hsap.sh
│   └── run_test_mmus.sh
├── docs/
│   ├── usage.md
│   ├── output.md
│   └── nf-core-reorganization.md
├── assets/
│   └── multiqc_config.yaml
└── meta/
    └── mappings/               # Moved from meta/ root
        ├── census_map_human.tsv
        ├── census_map_mouse_author.tsv
        ├── cell_type_markers.tsv
        └── color_mapping.tsv
```

---

## Modules Created (12 total)

| Module | Description | Original Process |
|--------|-------------|------------------|
| `save_params` | Write run parameters to YAML | `save_params_to_file` |
| `run_setup` | Download SCVI model from Census | `runSetup` |
| `get_census_adata` | Fetch reference data from Census | `getCensusAdata` |
| `map_query` | Process queries through SCVI model | `mapQuery` |
| `query_process_seurat` | Prepare queries for Seurat | `queryProcessSeurat` |
| `ref_process_seurat` | Prepare references for Seurat | `refProcessSeurat` |
| `rf_predict` | SCVI Random Forest prediction | `rfPredict` |
| `predict_seurat` | Seurat label transfer | `predictSeurat` |
| `classify_all` | Compute metrics (F1, confusion) | `classifyAll` |
| `plot_f1_results` | Plot F1 score distributions | `plotF1ResultsAdata/Seurat` |
| `plot_qc_combined` | Generate QC plots | `plotQC_combined` |
| `multiqc` | Create MultiQC HTML report | `runMultiQC` |

---

## Subworkflows Created (4 total)

| Subworkflow | Description | Modules Used |
|-------------|-------------|--------------|
| `prepare_references` | Downloads SCVI model and Census references | `run_setup`, `get_census_adata`, `ref_process_seurat` |
| `scvi_pipeline` | SCVI-based prediction workflow | `rf_predict`, `classify_all` |
| `seurat_pipeline` | Seurat label transfer workflow | `query_process_seurat`, `predict_seurat`, `classify_all` |
| `qc_reporting` | QC visualization and reporting | `plot_qc_combined`, `multiqc` |

---

## Configuration Changes

### Split into Multiple Files

**Before:** Single `nextflow.config` with all settings

**After:**
- `conf/base.config` - Process defaults with resource labels (`process_low`, `process_medium`, `process_high`)
- `conf/params.config` - All default parameter values
- `conf/modules.config` - Per-module publishDir configurations
- `conf/test_*.config` - Test-specific parameter overrides

### New Profiles Added

| Profile | Purpose |
|---------|---------|
| `conda` | Enable conda environments |
| `slurm` | Execute on SLURM cluster |
| `test` | Generic subsampled test |
| `test_hsap` | Human-specific test (100 query cells) |
| `test_mmus` | Mouse-specific test (100 query cells) |
| `debug` | Enable debug logging |

### Pipeline Manifest Added

```groovy
manifest {
    name            = 'cell-annotation-benchmark'
    author          = 'rschwartz'
    description     = 'Benchmarking pipeline for cell type annotation methods'
    mainScript      = 'main.nf'
    nextflowVersion = '!>=23.04.0'
    version         = '1.0.0'
}
```

---

## Test Scripts

**Before:** `scripts/run_hsap_test.sh` and `scripts/run_mmus_test.sh` with inline parameters

**After:** Simplified test runners using profiles:
```bash
# tests/run_test_hsap.sh
nextflow run main.nf -profile conda,test_hsap -work-dir hsap -resume "$@"

# tests/run_test_mmus.sh
nextflow run main.nf -profile conda,test_mmus -work-dir mmus -resume "$@"
```

---

## Files Moved/Reorganized

| Original Location | New Location |
|-------------------|--------------|
| `meta/census_map_*.tsv` | `meta/mappings/` |
| `meta/cell_type_markers.tsv` | `meta/mappings/` |
| `meta/color_mapping.tsv` | `meta/mappings/` |
| `meta/multiqc_config.yaml` | `assets/` |
| `scripts/run_*_test.sh` | `tests/` (rewritten) |

---

## Files Deleted

- `scripts/run_hsap_test.sh` (replaced by `tests/run_test_hsap.sh`)
- `scripts/run_mmus_test.sh` (replaced by `tests/run_test_mmus.sh`)
- `.github/copilot-instructions.md`

---

## How to Run

```bash
# Basic run
nextflow run main.nf -profile conda

# With parameters file
nextflow run main.nf -profile conda -params-file params.json

# Run tests
./tests/run_test_hsap.sh
./tests/run_test_mmus.sh
```

---

## Commits

1. `f76f316` - Reorganize pipeline to nf-core best practices
2. `5156efe` - Move test scripts to tests/ with proper nf-core profiles
3. `8bffa6a` - Update README and fix modules.config warnings