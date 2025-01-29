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
import json

# Function to parse command line arguments
def parse_arguments():
  parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
  parser.add_argument('--model_path', type=str, default="/space/grp/rschwartz/rschwartz/biof501_proj/scvi-human-2024-07-01", help='Path to the scvi model file')
  parser.add_argument('--subsample_query', type=int, help='Number of cells to subsample from the query')
  parser.add_argument('--relabel_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/pineda_relabel.tsv.gz")
  parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/queries/pineda.h5ad")
  parser.add_argument('--batch_key', type=str, default="sample")
  parser.add_argument('--join_key', type=str, default="")
  #parser.add_argument('--split_sex', action='store_true', help='Split by sex')
  #parser.add_argument('--split_disease', action='store_true', help='Split by disease')
  parser.add_argument('--ref_keys', type=str, nargs='+', default=["rachel_subclass", "rachel_class", "rachel_family"])
  parser.add_argument('--remove_unknown', action='store_true', help='Remove cells with unknown labels')
  parser.add_argument('--seed', type=int, default=42)
   
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

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
  query = ad.read_h5ad(query_path)

  
  query = relabel(query,relabel_path=relabel_path, join_key=join_key,sep="\t")

  if subsample_query:
    query= query[np.random.choice(query.n_obs, size=subsample_query, replace=False), :] if query.n_obs > subsample_query else query
    query.obs.index = query.obs.index.astype(str)
 
  if args.remove_unknown:
    query = query[query.obs[ref_keys[0]] != "unknown"]
    
  query.obs.index = range(query.n_obs)
  
  query = process_query(query, model_path, batch_key)
  new_query_name = os.path.basename(query_path).replace(".h5ad","_processed")
  query.write_h5ad(f"{new_query_name}.h5ad")
  
  
  
if __name__ == "__main__":
    main()