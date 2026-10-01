#!/usr/bin/env python
"""Rebuild per-sample metrics at any probability cutoff from a cutoff-0 run, without rerunning the pipeline.

Run: /home/rschwartz/anaconda3/envs/scanpyenv/bin/python scripts/rebuild_cutoff_metrics.py --cutoffs 0.25
At cutoff 0 every cell keeps its top prediction (predicted_<key> in *.predictions.0.0.tsv.gz) and its confidence.
For a cutoff c, classify_raw calls a cell unknown when confidence <= c, and unknown propagates to every level
(aggregate_labels fills unmapped labels with the lower label; map_valid_labels copies the higher-level prediction),
so setting all predicted_<key> to "unknown" for those cells and rescoring with utils.evaluate_sample_predictions
gives what CLASSIFY_ALL would have written. Author-unlabeled cells are not in the predictions files (they are
dropped before scoring); their unknown rate comes from *.unlabeled_cells.0.0.tsv.gz the same way.
Writes rebuilt_metrics.tsv (one row per method, query, reference, key, cutoff), rebuilt_unlabeled_unknown.tsv,
rebuilt_per_label.tsv (subclass counts per query and cutoff: true cells, predicted cells, correct cells, true cells called unknown, and the
unlabeled cells whose top prediction is that label, with how many of them are unknown) and manifest.json to --outdir. With --validate <cutoff dir of a real run at that cutoff> it also writes
validation.tsv comparing the rebuilt weighted F1, accuracy, NMI and ARI with that run's summary.scores files.
"""
import argparse, glob, os, sys, json, hashlib, datetime, subprocess
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{ROOT}/bin")
BASE = "/cosmos/data/nextflow-eval-pipeline/results/census-map-fixes_with_unlabeled_test/mus_musculus/ratio_native/ref_500"
REFS = "/cosmos/data/nextflow-eval-pipeline/results/cache/refs/census-map-fixes/mus_musculus/2024-07-01/brain/dataset_id/sub_500/refs"
METHODS = ["scvi_knn", "scvi_rf", "seurat"]
KEYS = ["subclass", "class", "family", "global"]
SCORES = ["weighted_f1", "macro_f1", "overall_accuracy", "nmi", "ari"]

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--run0", default=f"{BASE}/cutoff_0", help="pipeline output dir of the cutoff 0 run")
ap.add_argument("--refs", default=REFS, help="dir with <reference>.ref_counts.tsv")
ap.add_argument("--mapping", default=f"{ROOT}/meta/census_map_mouse_author.tsv")
ap.add_argument("--cutoffs", type=float, nargs="+", default=[round(x, 2) for x in np.arange(0, 0.96, 0.05)])
ap.add_argument("--validate", default=None, help="output dir of a real run at one of --cutoffs")
ap.add_argument("--outdir", default=None, help="default: <run0>/rebuilt_cutoffs")
ap.add_argument("--workers", type=int, default=8)
ap.add_argument("--limit", type=int, default=None, help="only the first N predictions files (smoke test)")
args = ap.parse_args()
OUT = args.outdir or f"{args.run0}/rebuilt_cutoffs"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def one(task):
    """Rescore one (method, study, reference, query) at every cutoff."""
    from utils import evaluate_sample_predictions
    pred_path, cells_path, method, study, ref, query = task
    mapping = pd.read_csv(args.mapping, sep="\t")
    rc = pd.read_csv(f"{args.refs}/{ref}.ref_counts.tsv", sep="\t")
    lookup = {k: g.set_index("label")["ref_support"].to_dict() for k, g in rc.groupby("key")}
    pred = pd.read_csv(pred_path, sep="\t")
    cells = pd.read_csv(cells_path, sep="\t")
    cells = cells[cells["author_unlabeled"].astype(str) == "True"]
    rows, urows, lrows = [], [], []
    for c in args.cutoffs:
        q = pred.copy()
        unk = (q["confidence"] <= c).values
        for k in KEYS:
            q[f"predicted_{k}"] = q[f"predicted_{k}"].astype(object)
            q.loc[unk, f"predicted_{k}"] = "unknown"
        m = evaluate_sample_predictions(q, KEYS, mapping, ref_counts_lookup=lookup)
        for k in KEYS:
            r = {"method": method, "study": study, "reference": ref, "query": query, "key": k, "cutoff": c,
                 "n_labeled": len(q), "n_labeled_unknown": int(unk.sum()),
                 "weighted_f1": m[k]["weighted_metrics"].get("f1_score"), "macro_f1": m[k]["macro_metrics"].get("f1_score"),
                 "overall_accuracy": m[k]["overall_accuracy"], "nmi": m[k]["nmi"], "ari": m[k]["ari"]}
            rows.append(r)
        true, pr = q["subclass"].astype(str), q["predicted_subclass"].astype(str)
        n_true, n_pred, tp = true.value_counts(), pr.value_counts(), true[true == pr].value_counts()
        n_tunk = true[pr == "unknown"].value_counts()
        n_upred = cells["predicted_subclass"].astype(str).value_counts()
        n_uunk = cells.loc[cells["confidence"] <= c, "predicted_subclass"].astype(str).value_counts()
        for lab in sorted(set(n_true.index) | set(n_upred.index)):
            lrows.append({"method": method, "study": study, "reference": ref, "query": query, "cutoff": c, "label": lab,
                          "n_true": int(n_true.get(lab, 0)), "n_pred": int(n_pred.get(lab, 0)), "tp": int(tp.get(lab, 0)), "n_true_unknown": int(n_tunk.get(lab, 0)),
                          "n_unlabeled_pred": int(n_upred.get(lab, 0)), "n_unlabeled_pred_unknown": int(n_uunk.get(lab, 0))})
        urows.append({"method": method, "study": study, "reference": ref, "query": query, "cutoff": c,
                      "n_unlabeled": len(cells), "n_unlabeled_unknown": int((cells["confidence"] <= c).sum())})
    return rows, urows, lrows


tasks = []
for m in METHODS:
    for p in sorted(glob.glob(f"{args.run0}/{m}/*/*/*/predicted_meta/*.predictions.0.0.tsv.gz")):
        parts = os.path.relpath(p, args.run0).split("/")
        _, study, ref, query = parts[:4]
        cells = glob.glob(f"{os.path.dirname(os.path.dirname(p))}/unlabeled_qc/*.unlabeled_cells.0.0.tsv.gz")
        assert len(cells) == 1, (p, cells)
        tasks.append((p, cells[0], m, study, ref, query))
if args.limit:
    tasks = tasks[:args.limit]
print(len(tasks), "queries x", len(args.cutoffs), "cutoffs", flush=True)

rows, urows, lrows = [], [], []
with ProcessPoolExecutor(args.workers) as ex:
    for i, (r, u, l) in enumerate(ex.map(one, tasks)):
        rows += r; urows += u; lrows += l
        if i % 50 == 0: print(i, "done", flush=True)
os.makedirs(OUT, exist_ok=True)
res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/rebuilt_metrics.tsv", sep="\t", index=False)
pd.DataFrame(urows).to_csv(f"{OUT}/rebuilt_unlabeled_unknown.tsv", sep="\t", index=False)
pd.DataFrame(lrows).to_csv(f"{OUT}/rebuilt_per_label.tsv", sep="\t", index=False)

val = None
if args.validate:
    real = []
    for m in METHODS:
        for p in glob.glob(f"{args.validate}/{m}/*/*/*/label_transfer_metrics/*.summary.scores.tsv.gz"):
            d = pd.read_csv(p, sep="\t")
            d = d.groupby(["query", "reference", "key", "cutoff"], as_index=False)[SCORES].first()
            real.append(d.assign(method=m, study=os.path.relpath(p, args.validate).split("/")[1]))
    real = pd.concat(real)
    cut = float(real["cutoff"].iloc[0])
    both = res[res["cutoff"] == cut].merge(real, on=["method", "study", "reference", "key"], suffixes=("", "_real"))
    both = both[both["query"] == both["query_real"]]
    val = both[["method", "study", "reference", "query", "key"]].copy()
    for s in SCORES:
        val[f"{s}_rebuilt"], val[f"{s}_real"] = both[s].values, both[f"{s}_real"].values
        val[f"{s}_absdiff"] = (both[s] - both[f"{s}_real"]).abs().values
    val.to_csv(f"{OUT}/validation.tsv", sep="\t", index=False)
    print(f"validation at cutoff {cut}: {len(val)} rows of {len(res[res['cutoff'] == cut])} rebuilt")
    print(val[[f"{s}_absdiff" for s in SCORES]].describe().loc[["count", "max", "mean"]].round(6).to_string())

git = lambda *a: subprocess.run(["git", "-C", ROOT, *a], capture_output=True, text=True).stdout.strip()
json.dump({"timestamp": datetime.datetime.now().isoformat(timespec="seconds"), "script": os.path.abspath(__file__),
           "script_sha256": sha(os.path.abspath(__file__)), "git_branch": git("rev-parse", "--abbrev-ref", "HEAD"),
           "git_commit": git("rev-parse", "HEAD"), "git_dirty": bool(git("status", "--porcelain")),
           "data_version": {"run0": args.run0, "validate": args.validate, "mapping": args.mapping, "mapping_sha256": sha(args.mapping)},
           "cutoffs": args.cutoffs, "n_queries": len(tasks),
           "inputs_sha256": {t[0]: sha(t[0]) for t in tasks} | {t[1]: sha(t[1]) for t in tasks}},
          open(f"{OUT}/manifest.json", "w"), indent=1)
