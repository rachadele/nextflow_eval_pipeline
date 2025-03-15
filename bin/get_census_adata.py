#!/user/bin/python3

from pathlib import Path
import random
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
import yaml

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--subsample_ref', type=int, default=10)
    parser.add_argument('--relabel_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    parser.add_argument('--ref_collections', type=str, nargs = '+', default = [
        "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation",
        "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types",
        "Adult mouse cortical cell taxonomy revealed by single cell transcriptomics",
        "Tabula Muris Senis"
    ]) 
    parser.add_argument('--split_column', type=str, default="tissue")
    parser.add_argument('--ref_keys', type=str, nargs="+", default=["subclass","class","family"])
    parser.add_argument('--seed', type=int, default=42)
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
        
     
def create_ref_region_yaml(refs, outdir):
    ref_names = list(refs.keys())
    ref_region_yaml = {}
    for ref_name in ref_names:
        unique_regions = refs[ref_name].obs["tissue"].unique().tolist()
        
        # If there are multiple region types, handle them by including them as a list
        if len(unique_regions) == 1:
            ref_region_yaml[ref_name.replace(" ", "_").replace("\\/", "_") \
                .replace("(","").replace(")","") \
                .replace("\\", "") \
                .replace("'", "") \
                .replace(":", "") \
                .replace(";", "") \
                .replace("&", "") ] = unique_regions[0]
        else:
            ref_region_yaml[ref_name.replace(" ", "_").replace("\\/", "_") \
                .replace("(","").replace(")","") \
                .replace("\\", "") \
                .replace("'", "") \
                .replace(":", "") \
                .replace(";", "") \
                .replace("&", "") ] = "multiple regions"
    
    with open(os.path.join(outdir, "ref_region.yaml"), 'w') as file:
        yaml.dump(ref_region_yaml, file)
         
         
def main():
    # Parse command line arguments
    args = parse_arguments()

    # Set organism and census_version from arguments
    organism = args.organism
    census_version = args.census_version
    subsample_ref = args.subsample_ref
    relabel_path = args.relabel_path
    ref_collections = args.ref_collections
    split_column = args.split_column
    ref_keys = args.ref_keys
    SEED = args.seed

    refs = adata_functions.get_census(
        organism=organism,
        subsample=subsample_ref,
        split_column=split_column,
        census_version=census_version,
        relabel_path=relabel_path,
        ref_collections=ref_collections,
        seed=SEED,
        ref_keys=ref_keys,
    )
    
    if "All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2" in refs:
        refs.pop('All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2')
    
    
    print("finished fetching anndata")
    outdir = "refs"
    os.makedirs(outdir, exist_ok=True)

    for ref_name, ref in refs.items():
        if ref.shape[0] < 50:
            continue

        new_ref_name = (
            ref_name.replace(" ", "_")
            .replace("\\/", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("\\", "")
            .replace("'", "")
            .replace(":", "")
            .replace(";", "")
            .replace("&", "")
        )

        ref.write(os.path.join(outdir, f"{new_ref_name}.h5ad"))
        ref.obs.to_csv(os.path.join(outdir, f"{new_ref_name}.obs.tsv"), sep="\t")

    create_ref_region_yaml(refs, outdir)

      
if __name__ == "__main__":
    main()