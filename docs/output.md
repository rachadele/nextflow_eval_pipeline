# Output

## Directory structure

```
<outdir>/
├── refs/
│   ├── scvi/           # SCVI reference data
│   └── seurat/         # Seurat reference data (.rds)
├── scvi/
│   └── <study>/<ref>/<query>/
│       ├── label_transfer_metrics/
│       ├── confusion/
│       └── predicted_meta/
├── seurat/
│   └── <study>/<ref>/<query>/
│       ├── label_transfer_metrics/
│       ├── confusion/
│       └── predicted_meta/
├── multiqc_results/
│   └── multiqc_report.html
├── params.yaml
└── trace.txt
```

## Key outputs

- **label_transfer_metrics/**: F1 scores and classification metrics
- **confusion/**: Confusion matrices per label hierarchy
- **predicted_meta/**: Predicted cell type annotations
- **multiqc_report.html**: Interactive QC report
