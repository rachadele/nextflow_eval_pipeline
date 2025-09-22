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
import scvi
from sklearn.ensemble import RandomForestClassifier
import utils
from utils import *
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
from types import SimpleNamespace

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/cpsc545_proj/mapped_queries/velmeshev/whole_cortex/subsample_1000/query_mapped.h5ad")
    parser.add_argument('--ref_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/refs/whole_cortex.h5ad") #nargs ="+")
    parser.add_argument('--ref_keys', type=str, nargs='+', default=["subclass", "class", "family"])

    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

    
    
def main():
    SEED = 42
    random.seed(SEED)         # For `random`
    np.random.seed(SEED)      # For `numpy`
    # For `torch`'
    scvi.settings.seed = SEED # For `scvi`
    # Parse command line arguments
    args = parse_arguments()

    # Set variables from arguments
   # organism = args.organism
   # census_version = args.census_version
    query_path = args.query_path
    ref_path = args.ref_path
    ref_keys = args.ref_keys

    # Load query and reference datasets
    query = ad.read_h5ad(query_path)
    #query=process_query(query, model_path="")
    query_name = os.path.basename(query_path).replace("_processed.h5ad", "")
    ref = ad.read_h5ad(ref_path, backed="r")
    ref_name = os.path.basename(ref_path).replace(".h5ad", "")

    # Run classification and ROC analysis
    probs = rfc_pred(ref=ref, query=query, ref_keys=ref_keys, seed=SEED)
    
    # eventually make this a data frame and save to disk, then make ROC another script
    
    probabilities = probs[ref_keys[0]]['probabilities']
    class_labels = probs[ref_keys[0]]['class_labels']

    # Create a DataFrame
    prob_df = pd.DataFrame(probabilities, columns=class_labels)
    #save data frame to inteim probs/ dir
    outdir="probs"
    os.makedirs(outdir, exist_ok=True)
    prob_df.to_csv(os.path.join(outdir,f"{query_name}_{ref_name}.prob.df.tsv"),sep="\t", index=False,)
    query.obs.to_csv(f"{query_name}.obs.relabel.tsv", index=False, sep="\t")

if __name__ == "__main__":
    main()
    
