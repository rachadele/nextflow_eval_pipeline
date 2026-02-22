# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Pipeline Does

Benchmarks automated cell type annotation on scRNA-seq data using CellXGene Census references. It compares two strategies:
- **SCVI + Random Forest**: Embeds cells using a pre-trained scVI model, then trains RF classifiers on embeddings
- **Seurat Label Transfer**: Gaussian kernel on dual PCA projection of reference and query

Supports two organisms (`homo_sapiens`, `mus_musculus`) and multiple cell type granularity levels (`subclass`, `class`, `family`, `global`).

## Running the Pipeline

```bash
# Basic run
nextflow run main.nf -profile conda

# With a params file
nextflow run main.nf -profile conda -params-file params.mm.json

# Test profiles
nextflow run main.nf -profile conda,test_hsap
nextflow run main.nf -profile conda,test_mmus

# Or via test scripts
./tests/run_test_hsap.sh
./tests/run_test_mmus.sh
```

Example full run (mouse, SCTransform normalization):
```bash
nextflow run main.nf -params-file params.mm.json \
    --subsample_query 100 \
    --subsample_ref 500 \
    --ref_split dataset_id \
    --cutoff 0.05 \
    --normalization_method SCT \
    --use_gap false \
    --outdir_prefix "results/mmus/100/dataset_id/SCT"
```

Work directory is set to `/cosmos/data/nextflow-eval-pipeline-work` in `nextflow.config`.

## Key Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `--organism` | null | `homo_sapiens` or `mus_musculus` |
| `--census_version` | null | CellXGene Census version (e.g. `2024-07-01`, `2025-01-30`) |
| `--subsample_ref` | 500 | Cells per cell type in reference |
| `--subsample_query` | 100 | Total cells in query (null = all) |
| `--ref_keys` | `["subclass","class","family","global"]` | Cell type hierarchy levels |
| `--cutoff` | 0 | Probability threshold for predictions |
| `--use_gap` | true | Gap-based thresholding |
| `--normalization_method` | `SCT` | Seurat normalization method |
| `--batch_correct` | true | Enable batch correction |
| `--ref_collections` | null | Reference dataset names from Census |
| `--queries_adata` | null | Glob path to query h5ad files |
| `--relabel_q` | null | Glob path to query relabeling files |
| `--relabel_r` | null | Reference relabeling mapping |

All defaults live in `conf/params.config`.

## Architecture

### Entry Point & Subworkflows

`main.nf` orchestrates four subworkflows in `subworkflows/local/`:
1. `prepare_references` — Downloads SCVI model + fetches Census reference data
2. `scvi_pipeline` — SCVI embedding → RF prediction → classification metrics
3. `seurat_pipeline` — Seurat preprocessing → label transfer → classification metrics
4. `qc_reporting` — QC plots + MultiQC HTML report

### Modules

Each process is a standalone module in `modules/local/<name>/main.nf`. The 12 modules map to scripts in `bin/`:

| Module | Script | Purpose |
|--------|--------|---------|
| `run_setup` | `setup.py` | Download SCVI model from Census |
| `get_census_adata` | `get_census_adata.py` | Fetch & process Census reference h5ad |
| `map_query` | `process_query.py` | Map query through SCVI model |
| `rf_predict` | `predict_scvi.py` | RF prediction on SCVI embeddings |
| `ref_process_seurat` | `ref_preprocessing.R` | Seurat reference preprocessing |
| `query_process_seurat` | `seurat_preprocessing.R` | Seurat query preprocessing |
| `predict_seurat` | `predict_seurat.R` | Seurat label transfer |
| `classify_all` | `classify_all.py` | F1 scores + confusion matrices |
| `plot_qc_combined` | `plot_QC_combined.py` | QC visualizations |
| `multiqc` | MultiQC | Aggregate QC HTML report |
| `plot_f1_results` | — | Plot F1 distributions |
| `save_params` | — | Save run parameters to YAML |

### Shared Utilities

- `bin/utils.py` (~44 KB) — shared Python utilities (cell type mapping, relabeling, metrics)
- `bin/seurat_functions.R` (~42 KB) — shared R/Seurat utilities

### Metadata

`meta/` contains cell type hierarchy mappings and relabeling files:
- `census_map_human.tsv` / `census_map_mouse_author.tsv` — cell type hierarchy mappings
- `relabel_homo_sapiens/` / `relabel_mus_musculus/` — per-study relabeling files that harmonize author labels to a common hierarchy

### Conda Environments

- Python (scanpy, scvi-tools, scikit-learn): `/home/rschwartz/anaconda3/envs/scanpyenv`
- R 4.3 (Seurat): `/home/rschwartz/anaconda3/envs/r4.3`

### Resource Labels (SLURM)

Defined in `conf/base.config`:
- `process_low`: 4 GB, 1 h
- `process_medium`: 36 GB, 8 h
- `process_high`: 72 GB, 16 h
- `process_long`: 20 h (time only)

SLURM constraint: `-C thrd64 --cpus-per-task=20`

## Output Structure

```
<outdir>/
├── refs/scvi/          # SCVI model + reference h5ads + UMAP plots
├── refs/seurat/        # Seurat RDS reference objects
├── scvi/<study>/<ref>/<query>/
│   ├── label_transfer_metrics/   # F1 scores (TSV)
│   ├── confusion/                # Confusion matrices
│   └── predicted_meta/           # Predicted cell type metadata
├── seurat/<study>/<ref>/<query>/
│   └── (same structure as scvi/)
├── multiqc_results/multiqc_report.html
└── params.yaml
```
