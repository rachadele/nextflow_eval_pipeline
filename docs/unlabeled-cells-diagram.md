# Unlabeled-cell cutoff test: logic

Goal: test whether a classifier probability cutoff calls author-unlabeled cells "unknown" without hurting labeled cells.
Style and colors follow `worfklow-diagrams/workflow_diagram.mmd` (blue = data build, green = benchmark, amber = scoring).

```mermaid
flowchart TD
    IN_MEX(["Gemma mex + author CTA\n(per sample)"])
    IN_CENSUS(["Census reference caches\n(git_branch = census-map-fixes)"])

    subgraph B["① build_unlabeled_queries.py  ·  Query Construction"]
        direction TB
        AUG["Add every barcode missing from the CTA\ncell_type = author_unlabeled"]
        REGEN["Regenerate h5ads\n(regenerate_h5ads_from_local_cta.sh)"]
        H5["Query h5ads\nobs bool author_unlabeled\nnative unlabeled:labeled ratio"]
        AUG --> REGEN --> H5
    end

    subgraph P["② nextflow_eval_pipeline  ·  Benchmark at cutoff 0.25"]
        direction TB
        SUB["process_query.py\nsubsample labeled cells to 100\nkeep all unlabeled cells of the sample"]
        QC["Scrublet + MAD QC flags\nQC plots (unlabeled shown as unscored)"]
        MAP["Map queries through scVI model"]
        CLF["Classify\nscVI + random forest, Seurat label transfer"]
        SUB --> QC
        SUB --> MAP --> CLF
    end

    subgraph S["③ classify_all.py  ·  Scoring"]
        direction TB
        SPLIT{"author_unlabeled?"}
        MET["Labeled cells only\nF1, confusion, NMI, ARI"]
        UNK["Unlabeled cells only\nper-cell table + unknown rate\n(unlabeled_qc/)"]
        CMP["Unknown rate of unlabeled cells\nvs F1 cost on labeled cells"]
        SPLIT -- no --> MET --> CMP
        SPLIT -- yes --> UNK --> CMP
    end

    OUT(["Cutoff sweep summary\nand QC plots"])

    IN_MEX --> AUG
    IN_CENSUS --> MAP
    H5 -- "query h5ads" --> SUB
    CLF -- "predictions + probabilities" --> SPLIT
    QC --> OUT
    CMP --> OUT

    classDef bnode fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    classDef pnode fill:#dcfce7,stroke:#22c55e,color:#14532d
    classDef snode fill:#fef3c7,stroke:#f59e0b,color:#78350f
    classDef ioNode fill:#f8fafc,stroke:#64748b,stroke-dasharray:5 5,color:#1e293b

    class AUG,REGEN,H5 bnode
    class SUB,QC,MAP,CLF pnode
    class SPLIT,MET,UNK,CMP snode
    class IN_MEX,IN_CENSUS,OUT ioNode
```

Unlabeled cells never enter F1, confusion, NMI or ARI. They are scored only by the fraction called "unknown" at the cutoff.
