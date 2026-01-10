# Cell Type Annotation Benchmarking Pipeline

A Nextflow DSL2 pipeline for benchmarking automated cell type annotation methods on scRNA-seq data using CellXGene Census references.

## Table of Contents
- [Description](#description)
- [Pipeline Structure](#pipeline-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Output](#output)
- [Methods](#methods)
- [References](#references)

---

## Description

Single-cell RNA sequencing (scRNA-seq) provides crucial insights into cell-type-specific gene expression, particularly in the context of disease. However, meta-analysis of scRNA-seq datasets remains challenging due to inconsistent and often absent cell-type annotations across publicly available repositories, such as the Gene Expression Omnibus (GEO).

This pipeline benchmarks two prominent cell type annotation strategies:
- **SCVI + Random Forest**: A Random Forest classifier trained on cell embeddings using single-cell variational inference (scVI)
- **Seurat Label Transfer**: A Gaussian kernel trained on PCA projection of query onto reference datasets

Reference data comprises 2 studies, 10 brain regions, 12 individual dissections, and 3 levels of cell type granularity from CellXGene Census.

---

## Pipeline Structure

```
nextflow_eval_pipeline/
├── main.nf                 # Main workflow orchestrator
├── nextflow.config         # Pipeline configuration
├── conf/
│   ├── base.config         # Process defaults and resource labels
│   ├── params.config       # Default parameters
│   ├── modules.config      # Per-module publish settings
│   ├── test_hsap.config    # Human test configuration
│   └── test_mmus.config    # Mouse test configuration
├── modules/local/          # Individual process definitions
│   ├── run_setup/          # Download SCVI model
│   ├── get_census_adata/   # Fetch Census reference data
│   ├── map_query/          # Process queries through SCVI
│   ├── rf_predict/         # SCVI Random Forest prediction
│   ├── predict_seurat/     # Seurat label transfer
│   ├── classify_all/       # Compute metrics
│   └── ...
├── subworkflows/local/     # Logical groupings
│   ├── prepare_references/ # Download model + Census refs
│   ├── scvi_pipeline/      # SCVI prediction workflow
│   ├── seurat_pipeline/    # Seurat prediction workflow
│   └── qc_reporting/       # QC plots + MultiQC
├── bin/                    # Python and R scripts
├── meta/
│   ├── mappings/           # Census maps, markers, colors
│   └── relabel_*/          # Dataset relabeling files
├── tests/                  # Test runner scripts
├── docs/                   # Documentation
└── assets/                 # Static assets (MultiQC config)
```

---

## Installation

### Requirements
- Nextflow >= 23.04.0
- Conda (for environment management)
- Python 3.8+ with scanpy, scvi-tools
- R 4.3+ with Seurat

### Setup
```bash
git clone https://github.com/rachadele/nextflow_eval_pipeline.git
cd nextflow_eval_pipeline
```

---

## Usage

### Basic Run
```bash
nextflow run main.nf -profile conda
```

### With Parameters File
```bash
nextflow run main.nf -profile conda -params-file params.json
```

### Run Tests
```bash
# Human test (subsampled)
./tests/run_test_hsap.sh

# Mouse test (subsampled)
./tests/run_test_mmus.sh

# Or using profiles directly
nextflow run main.nf -profile conda,test_hsap
nextflow run main.nf -profile conda,test_mmus
```

### Available Profiles

| Profile | Description |
|---------|-------------|
| `conda` | Use conda environments |
| `slurm` | Execute on SLURM cluster |
| `test` | Generic subsampled test |
| `test_hsap` | Human-specific test (100 cells) |
| `test_mmus` | Mouse-specific test (100 cells) |
| `debug` | Enable debug logging |

---

## Configuration

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--organism` | Species (`homo_sapiens` or `mus_musculus`) | `homo_sapiens` |
| `--census_version` | CellXGene Census version | `2025-01-30` |
| `--ref_collections` | Reference collections from Census | See config |
| `--ref_keys` | Label hierarchy levels | `["subclass", "class", "family", "global"]` |
| `--subsample_ref` | Cells per type in reference | `500` |
| `--subsample_query` | Total cells in query (null = all) | `null` |
| `--normalization_method` | Seurat normalization | `SCT` |
| `--cutoff` | Probability threshold | `0` |
| `--use_gap` | Use gap-based thresholding | `true` |
| `--outdir` | Output directory | Computed from params |

### Configuration Files

- `conf/params.config` - Edit default parameters
- `conf/base.config` - Modify resource allocations
- `conf/test_*.config` - Customize test runs

---

## Output

```
<outdir>/
├── refs/
│   ├── scvi/           # SCVI reference data (.h5ad)
│   └── seurat/         # Seurat reference data (.rds)
├── scvi/
│   └── <study>/<ref>/<query>/
│       ├── label_transfer_metrics/   # F1 scores
│       ├── confusion/                # Confusion matrices
│       └── predicted_meta/           # Predictions
├── seurat/
│   └── <study>/<ref>/<query>/
│       ├── label_transfer_metrics/
│       ├── confusion/
│       └── predicted_meta/
├── multiqc_results/
│   └── multiqc_report.html   # Interactive QC report
├── params.yaml               # Run parameters
└── trace.txt                 # Execution trace
```

---

## Methods

### SCVI + Random Forest
A Random Forest classifier trained on latent embeddings from a pre-trained single-cell variational autoencoder (scVI) model from CellXGene Census.

### Seurat Label Transfer
A Gaussian kernel trained on dual PCA projection of reference and query datasets using the Seurat package.

---

## References

- Lopez, R., et al. "Deep Generative Modeling for Single-Cell Transcriptomics." Nature Methods, 2018.
- Jorstad, N.L., et al. "Transcriptomic Cytoarchitecture Reveals Principles of Human Neocortex Organization." Science, 2023.
- Abdulla, S., et al. "CZ CELL×GENE Discover: A Single-Cell Data Platform." bioRxiv, 2023.
- Pasquini, G., et al. "Automated methods for cell type annotation on scRNA-seq data." Computational and Structural Biotechnology Journal, 2021.
- Lotfollahi, M., et al. "The Future of Rapid and Automated Single-Cell Data Analysis Using Reference Mapping." Cell, 2024.

---

![workflow DAG](dag.png)
