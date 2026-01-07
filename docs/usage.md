# Usage

## Basic execution

```bash
nextflow run main.nf -profile conda
```

## With parameters file

```bash
nextflow run main.nf -profile conda -params-file params.json
```

## Available profiles

- `conda` - Use conda environments
- `slurm` - Execute on SLURM cluster
- `test` - Run with subsampled data for testing
- `debug` - Enable debug mode with hostname logging

## Key parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--organism` | Species to analyze | `homo_sapiens` |
| `--census_version` | CellXGene Census version | `2025-01-30` |
| `--ref_collections` | Reference collections from Census | See config |
| `--subsample_ref` | Cells per type in reference | `500` |
| `--subsample_query` | Total cells in query (null = all) | `null` |
| `--normalization_method` | Seurat normalization | `SCT` |
| `--outdir` | Output directory | Computed from params |
