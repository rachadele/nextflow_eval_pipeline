#!/usr/bin/env python3
"""
Fetch reference cell-type support counts from CellXGene Census and produce:
  - Per-organism/key pivot TSVs (rows=labels, cols=refs, values=ref_support)
"""

import warnings
warnings.filterwarnings("ignore")

import os
import re
import sys
import argparse

import pandas as pd

import utils
from utils import get_original_celltypes


# ---------------------------------------------------------------------------
# Census configs
# ---------------------------------------------------------------------------

CONFIGS = {
    "homo_sapiens": {
        "census_version": "2024-07-01",
        "organ": "brain",
        "split_column": "dataset_id",
        "relabel_path": "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv",
        "ref_collections": [
            "Transcriptomic cytoarchitecture reveals principles of human neocortex organization",
            "SEA-AD: Seattle Alzheimer's Disease Brain Cell Atlas",
            "Molecular and cellular evolution of the primate dorsolateral prefrontal cortex",
        ],
        "author_annotations_path": None,
        "original_celltype_columns": None,
    },
    "mus_musculus": {
        "census_version": "2024-07-01",
        "organ": "brain",
        "split_column": "dataset_id",
        "relabel_path": "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse_author.tsv",
        "ref_collections": [
            "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation",
            "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types",
            "Tabula Muris Senis",
        ],
        "author_annotations_path": "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/2024-07-01",
        "original_celltype_columns": "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/2024-07-01/original_celltype_columns.tsv",
    },
}

REF_KEYS = ["subclass", "class", "family", "global"]
SUBSAMPLE_REF = 500
SEED = 42

SKIP_REFS = {
    "All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Microglia - 24 months old wild-type and Rag1-KO",
    "Single-cell of aged oligodendrocytes",
}


# ---------------------------------------------------------------------------
# Reference display name shortening
# ---------------------------------------------------------------------------

STRIP_PREFIXES = [
    r"^Dissection:\s*",
    r"^Whole Taxonomy\s*[-–]\s*",
    r"Seattle Alzheimer['']?s Disease (Brain Cell )?Atlas\s*(SEA-AD)?",
    r"SEA-AD:\s*",
]

SHORT_NAMES = {
    "Human Multiple Cortical Areas SMART-seq": "Human MTG SMART-seq",
    "Single-nucleus transcriptome data from the dlPFC": "snRNA dlPFC",
    "whole cortex": "Whole Cortex",
    "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types": "Mouse MOp atlas",
    "Single-cell RNA-seq for all cortical & hippocampal regions (SMART-Seq v4)": "Mouse ctx+hpc SMART-seq",
    "Single-cell RNA-seq for all cortical & hippocampal regions (10x)": "Mouse ctx+hpc 10x",
    "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation": "Mouse isocortex+hpc taxonomy",
}


def shorten_name(name):
    if name in SHORT_NAMES:
        return SHORT_NAMES[name]
    result = name
    for pat in STRIP_PREFIXES:
        result = re.sub(pat, "", result, flags=re.IGNORECASE).strip()
    return result.strip(" -–") or name


def ref_sort_key(short_name):
    dissection_keywords = ["dissection", "angular", "prefrontal", "auditory",
                           "somatosensory", "visual", "cingulate", "dfc", "ang",
                           "acc", "s1", "v1", "a1"]
    is_dissection = int(any(k in short_name.lower() for k in dissection_keywords))
    return (is_dissection, short_name.lower())


# ---------------------------------------------------------------------------
# Data fetching
# ---------------------------------------------------------------------------

def fetch_refs(organism, cfg):
    print(f"\n=== {organism} ===", flush=True)

    original_celltypes = None
    if cfg["author_annotations_path"] and cfg["original_celltype_columns"]:
        original_celltypes = get_original_celltypes(
            columns_file=cfg["original_celltype_columns"],
            author_annotations_path=cfg["author_annotations_path"],
        )

    refs = utils.get_census(
        organism=organism,
        organ=cfg["organ"],
        subsample=SUBSAMPLE_REF,
        split_column=cfg["split_column"],
        census_version=cfg["census_version"],
        relabel_path=cfg["relabel_path"],
        ref_collections=cfg["ref_collections"],
        seed=SEED,
        ref_keys=REF_KEYS,
        original_celltypes=original_celltypes,
    )

    for skip in SKIP_REFS:
        refs.pop(skip, None)

    records = []
    for ref_name, ref in refs.items():
        if ref.shape[0] < 50:
            print(f"  Skipping {ref_name} (< 50 cells)", flush=True)
            continue
        print(f"  {ref_name}: {ref.shape[0]} cells", flush=True)
        for key in REF_KEYS:
            if key not in ref.obs.columns:
                continue
            for label, count in ref.obs[key].value_counts().items():
                records.append({
                    "organism": organism,
                    "ref_name": ref_name,
                    "key": key,
                    "label": label,
                    "ref_support": int(count),
                })

    return records


# ---------------------------------------------------------------------------
# Pivot + plot
# ---------------------------------------------------------------------------

def build_pivot(df, organism, key):
    sub = df[(df["organism"] == organism) & (df["key"] == key)]
    pivot = sub.pivot_table(index="label", columns="ref_name", values="ref_support", fill_value=0)
    pivot.index = sorted(pivot.index, key=lambda x: x.lower())
    pivot.columns = [shorten_name(c) for c in pivot.columns]
    return pivot[sorted(pivot.columns, key=ref_sort_key)]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Fetch ref support and plot coverage heatmaps")
    parser.add_argument("--organisms", nargs="+", default=["homo_sapiens", "mus_musculus"])
    parser.add_argument("--output_dir", default="results/ref_coverage")
    parser.add_argument("--keys", default=None, help="Comma-separated keys (default: all)")
    args = parser.parse_args()

    keys = [k.strip() for k in args.keys.split(",")] if args.keys else REF_KEYS
    os.makedirs(args.output_dir, exist_ok=True)

    all_records = []
    for organism in args.organisms:
        if organism not in CONFIGS:
            print(f"Unknown organism: {organism}", file=sys.stderr)
            continue
        all_records.extend(fetch_refs(organism, CONFIGS[organism]))

    df = pd.DataFrame(all_records)

    for organism in args.organisms:
        for key in keys:
            sub = df[(df["organism"] == organism) & (df["key"] == key)]
            if sub.empty:
                continue
            pivot = build_pivot(df, organism, key)

            tsv_path = os.path.join(args.output_dir, f"{organism}_{key}_ref_support.tsv")
            pivot.to_csv(tsv_path, sep="\t")
            print(f"  Saved TSV: {tsv_path}")


if __name__ == "__main__":
    main()
