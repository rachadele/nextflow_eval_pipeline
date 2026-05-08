#!/usr/bin/env python3
"""
Fetch CellXGene Census observation metadata for reference collections.

Replaces: get_census_obs.py, get_unique_ref_cells.py, get_original_celltypes.py, get_ref_support.py
Uses only the Census obs API — no h5ad downloads, no gene-expression fetch.

Outputs (all optional, at least one required):
  --obs_out     : raw obs TSV
  --counts_out  : cell-type count TSV
  --support_dir : per-key ref-support pivot TSVs (one file per ref_key)
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import json
import random
import argparse
import numpy as np
import pandas as pd
import cellxgene_census

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import utils
from utils import get_original_celltypes, map_author_labels


def fetch_obs(census_version, organism, organ, ref_collections):
    """Fetch filtered obs from Census (obs only, no expression data)."""
    with cellxgene_census.open_soma(census_version=census_version) as census:
        dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()
        value_filter = (
            f"tissue_general == '{organ}' and "
            f"is_primary_data == True and "
            f"disease == 'normal'"
        )
        obs = cellxgene_census.get_obs(census, organism, value_filter=value_filter)

    obs = obs.merge(dataset_info, on="dataset_id", suffixes=(None, "_y"))
    obs.drop(columns=[c for c in obs.columns if c.endswith("_y")], inplace=True)
    obs = obs[obs["collection_name"].isin(ref_collections)].copy()
    print(f"Fetched {len(obs)} cells for {organism}", flush=True)
    return obs


def subsample_obs(obs, subsample, ref_key, seed=42):
    """Subsample obs to at most `subsample` cells per label in ref_key, matching pipeline behavior."""
    random.seed(seed)
    np.random.seed(seed)
    keep = []
    for label, grp in obs.groupby(ref_key):
        ids = grp.index.tolist()
        keep.extend(random.sample(ids, min(len(ids), subsample)))
    return obs.loc[keep]


def apply_relabeling(obs, relabel_path):
    """Join obs with hierarchy mapping TSV to add subclass/class/family/global columns."""
    mapping = pd.read_csv(relabel_path, sep="\t")
    join_key = mapping.columns[0]
    if join_key not in obs.columns:
        raise ValueError(f"Relabel join key '{join_key}' not found in obs columns")
    obs = obs.merge(mapping, on=join_key, how="left", suffixes=(None, "_y"))
    obs.drop(columns=[c for c in obs.columns if c.endswith("_y")], inplace=True)
    return obs


def main():
    parser = argparse.ArgumentParser(
        description="Fetch Census obs metadata for reference collections (no h5ad download)."
    )
    parser.add_argument("--params_file", default=None,
                        help="JSON params file (e.g. params.hs.json); values are overridden by explicit flags")
    parser.add_argument("--organism", default=None,
                        help="homo_sapiens or mus_musculus")
    parser.add_argument("--organ", default=None)
    parser.add_argument("--census_version", default=None)
    parser.add_argument("--ref_collections", nargs="+", default=None,
                        help="Collection names to include")
    parser.add_argument("--relabel_path", default=None,
                        help="Hierarchy mapping TSV (adds subclass/class/family/global)")
    parser.add_argument("--original_celltype_columns", default=None,
                        help="TSV mapping dataset_title → author cell type column name")
    parser.add_argument("--author_annotations_path", default=None,
                        help="Directory containing per-dataset .obs.tsv files")
    parser.add_argument("--ref_keys", nargs="+",
                        default=["subclass", "class", "family", "global"],
                        help="Hierarchy levels to include in support pivots")
    parser.add_argument("--split_column", default="dataset_title",
                        help="Column used to split refs in support pivot output")
    parser.add_argument("--subsample_ref", type=int, default=None,
                        help="Max cells per subclass label (matches pipeline subsampling)")
    # outputs
    parser.add_argument("--obs_out", default=None,
                        help="Path for raw obs TSV")
    parser.add_argument("--counts_out", default=None,
                        help="Path for cell-type count TSV")
    parser.add_argument("--support_dir", default=None,
                        help="Directory for per-key ref-support pivot TSVs")
    args = parser.parse_args()

    # Load params file first, then let explicit CLI flags override
    params = {}
    if args.params_file:
        with open(args.params_file) as f:
            params = json.load(f)

    organism = args.organism or params.get("organism", "homo_sapiens")
    organ = args.organ or params.get("organ", "brain")
    census_version = args.census_version or params.get("census_version", "2025-01-30")
    ref_collections = args.ref_collections or params.get("ref_collections")
    relabel_path = args.relabel_path or params.get("relabel_r")
    original_celltype_columns = args.original_celltype_columns or params.get("original_celltype_columns")
    author_annotations_path = args.author_annotations_path or params.get("author_annotations_path")
    subsample_ref = args.subsample_ref or params.get("subsample_ref")

    if not ref_collections:
        parser.error("Provide --ref_collections or a --params_file with ref_collections")
    if not any([args.obs_out, args.counts_out, args.support_dir]):
        parser.error("Specify at least one output: --obs_out, --counts_out, or --support_dir")

    obs = fetch_obs(census_version, organism, organ, ref_collections)

    if original_celltype_columns and author_annotations_path:
        original_celltypes = get_original_celltypes(
            columns_file=original_celltype_columns,
            author_annotations_path=author_annotations_path,
        )
        obs = map_author_labels(obs, original_celltypes)

    if relabel_path:
        obs = apply_relabeling(obs, relabel_path)

    if args.obs_out:
        obs.to_csv(args.obs_out, sep="\t", index=False)
        print(f"Wrote {len(obs)} rows → {args.obs_out}")

    if args.counts_out:
        base_cols = ["cell_type", "dataset_title", "collection_name", "cell_type_ontology_term_id"]
        if "author_cell_type" in obs.columns:
            base_cols = ["author_cell_type"] + base_cols
        counts = obs[base_cols].value_counts().reset_index()
        counts.to_csv(args.counts_out, sep="\t", index=False)
        print(f"Wrote cell-type counts → {args.counts_out}")

    if args.support_dir:
        os.makedirs(args.support_dir, exist_ok=True)
        primary_key = args.ref_keys[0]  # subclass — subsampling is always by the finest level
        for key in args.ref_keys:
            if key not in obs.columns:
                print(f"  Skipping '{key}' (not in obs)", flush=True)
                continue
            records = []
            for ref_name, grp in obs.groupby(args.split_column):
                if subsample_ref and primary_key in grp.columns:
                    grp = subsample_obs(grp, subsample_ref, primary_key)
                for label, count in grp[key].value_counts().items():
                    records.append({"ref_name": ref_name, "label": label, "ref_support": int(count)})
            # Whole cortex: subsample from all obs combined
            whole = subsample_obs(obs, subsample_ref, primary_key) if subsample_ref and primary_key in obs.columns else obs
            for label, count in whole[key].value_counts().items():
                records.append({"ref_name": "whole cortex", "label": label, "ref_support": int(count)})
            if not records:
                continue
            pivot = (
                pd.DataFrame(records)
                .pivot_table(index="label", columns="ref_name", values="ref_support", fill_value=0)
                .sort_index()
            )
            pivot.columns = sorted(pivot.columns)
            path = os.path.join(args.support_dir, f"{organism}_{key}_ref_support.tsv")
            pivot.to_csv(path, sep="\t")
            print(f"  Wrote {path}")


if __name__ == "__main__":
    main()
