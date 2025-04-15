#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
import scvi
from scipy.sparse import csr_matrix
import warnings
import cellxgene_census
import cellxgene_census.experimental
from sklearn.ensemble import RandomForestClassifier
import adata_functions
from adata_functions import *
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import re
from scipy.stats import median_abs_deviation


# Function to parse command line arguments
def parse_arguments():
  parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
  parser.add_argument('--model_path', type=str, default="/space/grp/rschwartz/rschwartz/biof501_proj/scvi-human-2024-07-01", help='Path to the scvi model file')
  parser.add_argument('--subsample_query', type=int, help='Number of cells to subsample from the query')
  parser.add_argument('--relabel_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/pineda_relabel.tsv.gz")
  parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/queries/pineda.h5ad")
  parser.add_argument('--batch_key', type=str, default="sample")
  parser.add_argument('--join_key', type=str, default="")
  parser.add_argument('--ref_keys', type=str, nargs='+', default=["rachel_subclass", "rachel_class", "rachel_family"])
  parser.add_argument('--remove_unknown', action='store_true', help='Remove cells with unknown labels')
  parser.add_argument('--nmads', type=int, default=5, help='Number of median absolute deviations for outlier detection')
  parser.add_argument('--gene_mapping', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv", help='Path to the gene mapping file')
  parser.add_argument('--seed', type=int, default=42)
   
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

def is_outlier(adata, metric: str, nmads=5):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
  
def map_genes(query, gene_mapping):
    # Check if the "feature_name" column exists in query.var
    if "feature_name" not in query.var.columns:
        # Merge gene_mapping with query.var based on the index
        query.var = query.var.merge(gene_mapping["OFFICIAL_SYMBOL"], left_index=True, right_index=True, how="left")
        # Rename the merged column to "feature_name"
        query.var.rename(columns={"OFFICIAL_SYMBOL": "feature_name"}, inplace=True)
    return query


def get_qc_metrics(query, nmads):
    query.var["mito"] = query.var["feature_name"].str.startswith(("MT", "mt", "Mt"))
    query.var["ribo"] = query.var["feature_name"].str.startswith(("RP", "Rp", "rp"))
    query.var["hb"] = query.var["feature_name"].str.startswith(("HB", "Hb","hb"))
    # fill NaN values with False
    query.var["mito"].fillna(False, inplace=True)
    query.var["ribo"].fillna(False, inplace=True)
    query.var["hb"].fillna(False, inplace=True) 

    sc.pp.calculate_qc_metrics(query, qc_vars=["mito", "ribo", "hb"], log1p=True, inplace=True, percent_top=[20], use_raw=False)

    metrics = {
        "log1p_total_counts": "outlier_total_counts",
        "log1p_n_genes_by_counts": "outlier_n_genes_by_counts",
        "pct_counts_mito": "outlier_mito",
        "pct_counts_ribo": "outlier_ribo",
        "pct_counts_hb": "outlier_hb"
        }
    
    for metric, col_name in metrics.items():
        query.obs[col_name] = is_outlier(query, metric, nmads)
    
    query.obs["counts_outlier"] = query.obs["outlier_total_counts"] | query.obs["outlier_n_genes_by_counts"]
 
    return query



def main():

  # Parse command line arguments
  args = parse_arguments()
  SEED = args.seed
  random.seed(SEED)         # For `random`
  np.random.seed(SEED)      # For `numpy`
  scvi.settings.seed = SEED # For `scvi`
  # Set organism and census_version from arguments
  model_path = args.model_path
  subsample_query = args.subsample_query
  query_path = args.query_path
  relabel_path = args.relabel_path
  batch_key = args.batch_key
  join_key = args.join_key
  ref_keys = args.ref_keys
  nmads = args.nmads
  query = ad.read_h5ad(query_path)
  gene_mapping_path = args.gene_mapping
  gene_mapping = pd.read_csv(gene_mapping_path, sep="\t", header=0)
  # Drop rows with missing values in the relevant columns
  gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])
  # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
  gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
  gene_mapping.set_index("ENSEMBL_ID", inplace=True)
  query = map_genes(query, gene_mapping)
  query = relabel(query,relabel_path=relabel_path, join_key=join_key,sep="\t")

  if subsample_query:
    query= query[np.random.choice(query.n_obs, size=subsample_query, replace=False), :] if query.n_obs > subsample_query else query
    query.obs.index = query.obs.index.astype(str)
 
  if args.remove_unknown:
    query = query[query.obs[ref_keys[0]] != "unknown"]
    
  query.obs.index = range(query.n_obs)
  sc.pp.scrublet(query, batch_key=None)
  query = process_query(query, model_path, batch_key, seed=SEED)
  query = get_qc_metrics(query, nmads)
  new_query_name = os.path.basename(query_path).replace(".h5ad","_processed")
  query.write_h5ad(f"{new_query_name}.h5ad")
  
  
  
if __name__ == "__main__":
    main()