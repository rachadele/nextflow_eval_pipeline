# Unlabeled-cell cutoff test: logic

Goal: test whether a classifier probability cutoff labels author-unlabeled cells "unknown" without hurting labeled cells.

```mermaid
flowchart TD
    A[Gemma mex + author CTA per sample] --> B[build_unlabeled_queries.py]
    B -->|every barcode missing from the CTA gets cell_type = author_unlabeled| C[Augmented CTA + symlinked mex]
    C --> D[regenerate_h5ads_from_local_cta.sh]
    D --> E[Query h5ads with obs bool author_unlabeled<br/>native unlabeled:labeled ratio per sample]
    E --> F[process_query.py]
    F -->|subsample labeled cells to 100 as before<br/>add all unlabeled cells of the sample<br/>ground truth = author_unlabeled| G[Query Seurat objects]
    G --> H[QC: scrublet, MAD flags, plot_QC_combined.py]
    H -->|unlabeled cells plotted, predicted labels set to unscored| I[QC plots]
    G --> J[map_query / classifier, cutoff 0.25]
    R[Census reference caches<br/>git_branch = census-map-fixes] --> J
    J --> K[classify_all.py]
    K --> L{author_unlabeled?}
    L -->|no| M[F1, confusion, NMI/ARI metrics]
    L -->|yes| N[unlabeled_qc/ per-cell table<br/>unlabeled_unknown_rate]
    N --> O[Compare unknown rate of unlabeled cells<br/>vs F1 loss on labeled cells across cutoffs]
    M --> O
```

Unlabeled cells never enter F1, confusion, NMI or ARI. They are scored only by the fraction called "unknown" at the cutoff.
