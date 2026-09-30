# Handoff: add author-unlabeled cells to the mouse benchmark to test whether a cutoff catches them

From: `mask-annotation-overlap` (2026-09-29). Source notes: `/space/grp/rschwartz/rschwartz/mask-annotation-overlap/docs/UNLABELED_CELLS_QC.md` and `/space/grp/rschwartz/rschwartz/mask-annotation-overlap/README.md`.

## Background

The mouse query h5ads (`params.mm.json` → `queries_adata`, i.e. `get_gemma_data.nf/study_names_mouse.txt_author_true_process_samples_true/h5ad/`) hold only the cells the authors labeled. Gemma holds many more cells per sample. For example, GSE247339.2 has 140,868 cells in Gemma and 44,752 author labels.

The unlabeled cells are mostly low quality. In six of the seven studies, their median UMIs and detected genes fall well below those of labeled cells. In GSE247339.1/.2 the median unlabeled cell has under 130 UMIs. So the authors dropped these barcodes in QC before labeling, and the labels were not lost. Per-sample and per-study numbers are in `/space/grp/rschwartz/rschwartz/mask-annotation-overlap/results/qc_labeled_vs_unlabeled/`.

The sc-annotation pipeline's QC mask (per-sample MAD, nmads 20/5/5/5) flags only 14% of unlabeled cells in GSE124952 (1,133 of 7,956), from `/space/grp/rschwartz/rschwartz/mask-annotation-overlap/results/unannotated/unannotated_masked_counts.tsv`.

## Idea

Put some unlabeled cells back into the benchmark queries with a sentinel ground truth, then ask whether imposing a cutoff sends them to "unknown" more often than labeled cells. The cutoff can be the classifier probability `cutoff` (and `use_gap`), a QC floor, or both. Unlabeled cells should get flagged, and labeled cells should keep their labels.

## Proposed design

1. **Build augmented query h5ads in a separate dir** so the current queries and results stay untouched. For each study, write a CTA file that adds a row per unlabeled barcode (in the mex but absent from the author file) with `cell_type = "author_unlabeled"`. Then run `get_gemma_data.nf/bin/regenerate_h5ads_from_local_cta.sh` against a copy of the outdir holding those CTA files, for example `get_gemma_data.nf/study_names_mouse.txt_author_true_process_samples_true_with_unlabeled/`. Add an obs column `author_unlabeled` (bool) so the downstream code can split cells without relying on the label string.
2. **Choose which unlabeled cells to include.** Exclude GSE199460.2, where the author file covers only a vascular subset, so its unlabeled cells are not QC rejects. Also exclude samples with zero author labels: 2 in GSE185454, 1 in GSE199460.2 and 2 each in GSE247339.1/.2. Either include all unlabeled cells in the remaining samples or subsample them to a fixed ratio against labeled cells per sample. GSE247339 would otherwise be about 70% junk and could distort SCVI/SCT normalisation.
3. **Relabel tables.** Map `author_unlabeled` to itself at every level (subclass, class, family, global) in `meta/relabel_mus_musculus/*_relabel.tsv`, or have `classify_all.py` pass it through.
4. **Metrics.** Keep `author_unlabeled` cells out of the F1/precision/recall calculations, so existing metrics don't change. Add a separate output per study, sample, method, ref and cutoff: the fraction of unlabeled cells predicted "unknown" vs the fraction of labeled cells predicted "unknown". Also record per-cell UMIs and genes, so a QC floor can be scored on the same cells.
5. **Run.** Sweep `cutoff` over the usual values (0 up to 0.75) for scvi-knn/rf and Seurat on the mouse refs. Work on the `census-map-fixes` worktree (`nextflow_eval_pipeline-worktrees/census-map-fixes`, see `run_mmus_sweep_census_fixes.sh` for the current mouse sweep). Write to a new outdir so the existing census-fix sweep results are not overwritten.

## Open questions to settle before running

- Does any step (SCT, SCVI preprocessing, `process_query`) filter cells by minimum genes or UMIs? If so, many GSE247339 unlabeled cells drop out before classification, and that filter is itself a cutoff worth reporting.
- Probability cutoff, QC floor, or both? The user's goal is to mask all unlabeled cells. Earlier they suggested tuning nmads, so report both.
- Include all unlabeled cells or subsample them? Pick the ratio before running and record it in the run manifest.

## Status

Nothing has been built or run yet. The six-study sc-annotation-pipeline QC run (screen `mask-overlap-six`) is producing masks for the other mouse studies, which will give the MAD-mask comparison for the same cells.

## Progress (2026-09-29)

Code changes are in place in the `census-map-fixes` worktree, uncommitted. Nothing has been built or run. `params.mm.json` and `run_mmus_sweep_census_fixes.sh` were left alone.

### Answers to the open questions

**Does any step drop cells by genes/UMIs?** No step before classification filters query cells by a minimum gene or UMI count.

- `get_gemma_data.nf/bin/process_query_samples.py:106` drops cells with no `cell_type` after the CTA join. That is why the current h5ads hold only labeled cells. Giving unlabeled barcodes a CTA row keeps them. Lines 91-95 and 172-175 send whole samples with under 50 cells to `small_samples/`. There is no per-cell filter.
- `bin/process_query.py` randomly subsamples to `subsample_query` cells per sample (100 in the sweep; line 69 on this branch, line 59 before the edit), then runs scrublet (line 72). Scrublet filters `min_genes=3` on an internal copy only (scanpy 1.10.4 `_scrublet/__init__.py:191-192`). Cells under 3 genes stay in, with NaN `predicted_doublet`. `--remove_unknown` (line 92) drops only ground truth `"unknown"`. `get_qc_metrics` (line 99, `utils.py:831`) adds per-sample MAD flags (`umi_outlier`, `genes_outlier`, `counts_outlier`, `total_outlier`, nmads=5) and drops nothing. `utils.process_query` (`utils.py:290`) subsets genes to the scVI model's genes (`prepare_query_anndata`, line 304). It drops no cells.
- Seurat: `bin/seurat_preprocessing.R:23-36` converts with sceasy and runs SCTransform + PCA. It does not filter. The `min.features = 200` filter in `seurat_functions.R:13` (`process_sample`) and the MAD filters (`filter_valid_cells`, line 265) never run on the query path.
- `bin/classify_all.py` drops nothing. `utils.classify_raw` (line 585) and `classify_by_gap` (line 606) turn low-confidence cells into `"unknown"` one cell at a time.
- So the probability cutoff is the only thing that sends cells to "unknown". The MAD flags are stored and can serve as a QC floor offline. One risk: a near-empty GSE247339 barcode can have zero counts in the scVI model genes after line 304, and SCTransform may fail on a cell with zero UMIs. If QUERY_PROCESS_SEURAT fails on those samples, lower the ratio or drop zero-count cells in the build script.

**How does a ground truth missing from the relabel tables behave?** `utils.relabel` (`utils.py:97-113`) left-joins on `cell_type`. An unmapped label gets NaN `subclass`, and `process_query.py:83-89` raises a ValueError. So without a change, `author_unlabeled` would crash MAP_QUERY. `aggregate_labels` (`utils.py:117-133`) fills higher levels from the lower label when the census map lacks it (line 132), so a label set at subclass passes through class/family/global unchanged. `map_valid_labels` skips labels absent from the census map (`utils.py:326`). `classify_all.py` never checks ground-truth labels against the relabel tables.

**Where are metrics computed?** `utils.evaluate_sample_predictions` (`utils.py:671-751`), called from `classify_all.py:137`. It computes per-label precision/recall/F1 (line 703) and weighted/macro/micro averages (lines 733-741) over labels with nonzero query support and ref support. `classify_all.py` writes `label_transfer_metrics/*.summary.scores.tsv.gz`, which `evaluation_summary.nf` reads. If `author_unlabeled` stayed in, it would add its own label row and count as false negatives in micro averages, accuracy, NMI and ARI. It would also add to `total_cell_count`. So those cells are now removed before `map_valid_labels`/`evaluate_sample_predictions`.

**Probability cutoff vs QC floor:** report both. Each per-cell output row holds the prediction, confidence, UMIs, genes and MAD flags, so one sweep scores both.

**Ratio:** `--ratio` in the build script is required and has no default. Choose it before building. It goes into `manifest.json`.

### Files changed (worktree)

- `bin/process_query.py`: sets ground truth `author_unlabeled` at subclass after relabeling, so no relabel-table edits are needed (params.mm.json's `relabel_q` points at the main repo's tables anyway). When the h5ad has an `author_unlabeled` column, it draws the labeled subsample with the same RNG call as before and then adds unlabeled cells at the sample's ratio. The labeled cells then match runs without unlabeled cells. scVI latents and the RF/kNN predictions depend only on each cell, so labeled-cell metrics for scvi_rf/scvi_knn should match the census-fix sweep. Seurat metrics will shift, because SCTransform and anchors use the whole query.
- `bin/classify_all.py`: adds `--method` to the parser. It already came from the module but was ignored. When `author_unlabeled` is present, `write_unlabeled_unknown_rates()` writes `unlabeled_qc/<query>_<ref>.unlabeled_cells.<cutoff>.tsv.gz` (per cell: labels, prediction, confidence, `total_counts`, `n_genes_by_counts`, Seurat `nCount_RNA`/`nFeature_RNA`, mito %, MAD flags, `predicted_unknown`). It also writes `unlabeled_qc/<query>_<ref>.unlabeled_unknown_rate.<cutoff>.tsv.gz` (per group, labeled vs author_unlabeled: n, n_unknown, frac_unknown, frac_mad_outlier, median UMIs/genes, with study/query/ref/method/cutoff/use_gap). Those cells are then dropped from `predicted_meta`, confusion matrices and summary scores. Queries without the column behave exactly as before.
- `modules/local/classify_all/main.nf`: adds the optional `unlabeled_qc/**tsv.gz` output, published by the existing `**tsv.gz` pattern to `<outdir>/<method>/<study>/<ref>/<query>/unlabeled_qc/`.
- `scripts/build_unlabeled_queries.py` (new): builds the augmented CTAs and h5ads (see below).
- `scripts/run_mmus_sweep_with_unlabeled.sh` (new): copy of the census-fix sweep that sets `--queries_adata` and a new `--outdir`.

Known side effect: QC_REPORTING left-joins predictions onto the raw h5ad, so unlabeled cells show up in the QC plots with NaN predictions.

### Build script behaviour

`scripts/build_unlabeled_queries.py` only reads the source dir. It uses only samples that already have an h5ad there, which drops the zero-label and under-50 samples. It also skips any sample with zero author labels and refuses GSE199460.2. Per sample it draws `min(n_unlabeled, round(ratio * n_labeled))` unlabeled barcodes (seed 42). It writes `<study>.celltypes.tsv` with the extra rows and an `author_unlabeled` column, and symlinks the source mex sample dirs and `metadata_standardized/<study>` into the new dir. It then calls `regenerate_h5ads_from_local_cta.sh`, casts `obs["author_unlabeled"]` to bool, and checks that each h5ad's labeled count equals the source h5ad's `n_obs`. It writes `unlabeled_tally.tsv` and `manifest.json` (input sha256 and mtimes, output h5ad sha256, ratio, seed, git commit, timestamp). It refuses to run if the destination is not empty.

### Commands

(a) Build the augmented h5ads (pick the ratio first; 0.5 is only an example):

```bash
cd /space/grp/rschwartz/rschwartz/nextflow_eval_pipeline-worktrees/census-map-fixes
/home/rschwartz/anaconda3/envs/scanpyenv/bin/python scripts/build_unlabeled_queries.py --ratio 0.5
# output: /space/grp/rschwartz/rschwartz/get_gemma_data.nf/study_names_mouse.txt_author_true_process_samples_true_with_unlabeled/
```

(b) Launch the sweep into a new outdir (`/cosmos/data/nextflow-eval-pipeline/results/census-map-fixes_with_unlabeled/...`), inside screen/tmux:

```bash
cd /space/grp/rschwartz/rschwartz/nextflow_eval_pipeline-worktrees/census-map-fixes
bash scripts/run_mmus_sweep_with_unlabeled.sh 2>&1 | tee run_mmus_sweep_with_unlabeled.log
# confidence-gap variant:
USE_GAP=true bash scripts/run_mmus_sweep_with_unlabeled.sh 2>&1 | tee run_mmus_sweep_with_unlabeled_gap.log
```

Collect the results with `find /cosmos/data/nextflow-eval-pipeline/results/census-map-fixes_with_unlabeled -name '*.unlabeled_unknown_rate.*.tsv.gz'`.
