#!/usr/bin/env python3
"""Build mouse query h5ads that also hold author-unlabeled cells.

For each study, cells that are in the Gemma mex but absent from the author CTA
are added back with cell_type = "author_unlabeled", subsampled per sample to
--ratio unlabeled cells per labeled cell. The augmented CTA files go into a new
outdir that links to the source mex and sample metadata, and
get_gemma_data.nf/bin/regenerate_h5ads_from_local_cta.sh builds the h5ads there.
Each h5ad gets a bool obs column `author_unlabeled`.

Only samples that already have an h5ad in the source dir are used, so samples
with zero author labels (and those under 50 labeled cells) stay out. GSE199460.2
is always excluded: its author file covers a vascular subset only, so its
unlabeled cells are not QC rejects.

The source dir is only read. Run with the scanpyenv python:
    /home/rschwartz/anaconda3/envs/scanpyenv/bin/python scripts/build_unlabeled_queries.py --ratio 0.5
"""

import argparse
import datetime
import gzip
import hashlib
import json
import os
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd

GEMMA = "/space/grp/rschwartz/rschwartz/get_gemma_data.nf"
EXCLUDED_STUDIES = {"GSE199460.2"}
DEFAULT_STUDIES = ["GSE124952", "GSE181021.2", "GSE185454", "GSE214244.1", "GSE247339.1", "GSE247339.2"]


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--src", default=f"{GEMMA}/study_names_mouse.txt_author_true_process_samples_true")
    parser.add_argument("--dest", default=f"{GEMMA}/study_names_mouse.txt_author_true_process_samples_true_with_unlabeled")
    parser.add_argument("--ratio", type=float, required=True, help="unlabeled:labeled cells per sample (upper bound; all unlabeled cells if fewer)")
    parser.add_argument("--studies", nargs="+", default=DEFAULT_STUDIES)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--regenerate_script", default=f"{GEMMA}/bin/regenerate_h5ads_from_local_cta.sh")
    return parser.parse_args()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def mtime(path):
    return datetime.datetime.fromtimestamp(os.path.getmtime(path)).isoformat(timespec="seconds")


def read_barcodes(mex_dir):
    with gzip.open(os.path.join(mex_dir, "barcodes.tsv.gz"), "rt") as f:
        return [line.strip() for line in f]


def augment_study(study, src, dest, ratio, rng, inputs):
    cta_path = os.path.join(src, "cell_type_assignments", f"{study}.celltypes.tsv")
    cta = pd.read_csv(cta_path, sep="\t", dtype={"sample_id": str, "cell_id": str})
    inputs[cta_path] = {"sha256": sha256(cta_path), "mtime": mtime(cta_path)}

    h5ad_dir = os.path.join(src, "h5ad", study)
    samples = sorted(f[len(study) + 1:-len(".h5ad")] for f in os.listdir(h5ad_dir) if f.endswith(".h5ad"))

    os.makedirs(os.path.join(dest, "mex", study))
    new_rows, tally = [], []
    for sample_name in samples:
        sample_id = sample_name.split("_")[0]
        mex_dir = os.path.join(src, "mex", study, sample_name)
        barcodes = read_barcodes(mex_dir)
        inputs[os.path.join(mex_dir, "barcodes.tsv.gz")] = {"sha256": sha256(os.path.join(mex_dir, "barcodes.tsv.gz"))}
        labeled_ids = set(cta.loc[cta["sample_id"] == sample_id, "cell_id"])
        n_labeled = sum(b in labeled_ids for b in barcodes)
        unlabeled = [b for b in barcodes if b not in labeled_ids]
        if n_labeled == 0:
            print(f"  skip {study} {sample_name}: zero author labels")
            continue
        n_keep = min(len(unlabeled), int(round(ratio * n_labeled)))
        keep = sorted(rng.choice(len(unlabeled), size=n_keep, replace=False)) if n_keep else []
        new_rows += [(sample_id, unlabeled[i]) for i in keep]
        os.symlink(mex_dir, os.path.join(dest, "mex", study, sample_name))
        tally.append({"study": study, "sample": sample_name, "n_labeled": n_labeled,
                      "n_unlabeled_available": len(unlabeled), "n_unlabeled_included": n_keep})

    cta["author_unlabeled"] = False
    added = pd.DataFrame(new_rows, columns=["sample_id", "cell_id"])
    added["cell_type"] = "author_unlabeled"
    added["cell_type_uri"] = np.nan
    added["author_unlabeled"] = True
    out_cta = os.path.join(dest, "cell_type_assignments", f"{study}.celltypes.tsv")
    pd.concat([cta, added], ignore_index=True).to_csv(out_cta, sep="\t", index=False)

    meta_src = os.path.join(src, "metadata_standardized", study)
    os.symlink(meta_src, os.path.join(dest, "metadata_standardized", study))
    meta_file = os.path.join(meta_src, f"{study}_sample_meta_std.tsv")
    inputs[meta_file] = {"sha256": sha256(meta_file), "mtime": mtime(meta_file)}
    return tally


def main():
    args = parse_arguments()
    src, dest = os.path.abspath(args.src), os.path.abspath(args.dest)
    if src == dest:
        sys.exit("--dest must differ from --src")
    if os.path.exists(dest) and os.listdir(dest):
        sys.exit(f"{dest} exists and is not empty; remove it or pick another --dest")
    bad = EXCLUDED_STUDIES.intersection(args.studies)
    if bad:
        sys.exit(f"excluded studies requested: {sorted(bad)}")

    rng = np.random.default_rng(args.seed)
    for sub in ["mex", "cell_type_assignments", "metadata_standardized"]:
        os.makedirs(os.path.join(dest, sub), exist_ok=True)

    inputs, tally = {}, []
    for study in args.studies:
        print(f"{study}: writing augmented CTA")
        tally += augment_study(study, src, dest, args.ratio, rng, inputs)
    tally = pd.DataFrame(tally)

    subprocess.run(["bash", args.regenerate_script, dest, *args.studies], check=True)

    # regenerate writes every obs column as str; make author_unlabeled a bool and check the labeled counts
    outputs = {}
    tally["n_labeled_h5ad"], tally["n_unlabeled_h5ad"] = 0, 0
    for i, row in tally.iterrows():
        path = os.path.join(dest, "h5ad", row["study"], f"{row['study']}_{row['sample']}.h5ad")
        if not os.path.exists(path):
            print(f"  WARNING: no h5ad for {row['study']} {row['sample']}")
            continue
        adata = ad.read_h5ad(path)
        adata.obs["author_unlabeled"] = adata.obs["author_unlabeled"].astype(str) == "True"
        adata.write_h5ad(path)
        tally.loc[i, "n_unlabeled_h5ad"] = int(adata.obs["author_unlabeled"].sum())
        tally.loc[i, "n_labeled_h5ad"] = int((~adata.obs["author_unlabeled"]).sum())
        src_h5ad = os.path.join(src, "h5ad", row["study"], f"{row['study']}_{row['sample']}.h5ad")
        n_src = ad.read_h5ad(src_h5ad, backed="r").n_obs
        if n_src != tally.loc[i, "n_labeled_h5ad"]:
            print(f"  WARNING: {row['sample']} has {tally.loc[i, 'n_labeled_h5ad']} labeled cells, source h5ad has {n_src}")
        outputs[path] = sha256(path)
    tally.to_csv(os.path.join(dest, "unlabeled_tally.tsv"), sep="\t", index=False)

    git = lambda *a: subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(__file__)), *a],
                                    capture_output=True, text=True).stdout.strip()
    manifest = {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "script": os.path.abspath(__file__),
        "script_sha256": sha256(os.path.abspath(__file__)),
        "git_branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_commit": git("rev-parse", "HEAD"),
        "git_dirty": bool(git("status", "--porcelain")),
        "data_version": os.path.basename(src),
        "src": src,
        "dest": dest,
        "ratio_unlabeled_to_labeled": args.ratio,
        "seed": args.seed,
        "studies": args.studies,
        "excluded_studies": sorted(EXCLUDED_STUDIES),
        "regenerate_script": args.regenerate_script,
        "regenerate_script_sha256": sha256(args.regenerate_script),
        "inputs": inputs,
        "outputs_sha256": outputs,
        "totals": tally[["n_labeled_h5ad", "n_unlabeled_h5ad"]].sum().astype(int).to_dict(),
    }
    with open(os.path.join(dest, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(tally.groupby("study")[["n_labeled", "n_unlabeled_available", "n_unlabeled_included", "n_unlabeled_h5ad"]].sum())


if __name__ == "__main__":
    main()
